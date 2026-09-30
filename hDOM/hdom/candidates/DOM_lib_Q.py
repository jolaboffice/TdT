"""Shared primitives: the two correlation measures, TIFF-to-profile conversion, and small helpers.

Everything the pipeline measures reduces to two functions here, and they behave differently:

    corrcoef(x, y)   min-max rescales both inputs, then evaluates the Pearson formula. Pearson is
                     invariant to affine rescaling, so that first step changes nothing: the value
                     equals numpy.corrcoef(x, y)[0, 1] (agreement checked to ~1e-14). This is cc_r,
                     the red-channel score.
    dcc(x, y)        cosine of the angle between the raw vectors, WITHOUT subtracting the mean.
                     It is therefore NOT affine-invariant: adding an offset, or rescaling one
                     input, moves the value. This is cc_g, the green-channel score, and its
                     sensitivity to label brightness is a property of the measure, not a bug in
                     the caller.

Both are numba-compiled, so they are called once per candidate per shift without dominating runtime.

Profiles come from tif2RGprofile: each frame is collapsed to one value per column by taking the
column maximum, and the green profile is thresholded at mean + 0.1*(max - mean) with everything
below set to zero. That threshold is what "the labels pcc sees" means throughout.

Not called anywhere in the distributed pipeline, kept because older analysis scripts import them:
pcc, Dcc, crosscoef, low_pass_filter, compute_shift, cross_correlation_using_fft, two_images,
read_maps_file, show_image, show_image2, show_image3, show_image4.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from PIL import Image
import tifffile as tiff
from numpy.fft import fft, ifft, fftshift
import math
import numba

@numba.njit
def corrcoef(x,y):
    """Pearson correlation of x and y. cc_r.

    The min-max rescaling below is inherited from the original implementation and does not affect
    the result -- Pearson ignores affine transforms -- so this returns the same number as
    numpy.corrcoef. Returns 0 when either input is constant.
    """
    # normalize
    x_max=np.max(x)
    x_min=np.min(x)
    x_gap=x_max-x_min
    if (x_gap==0):
        x_gap=1
    x=(x-x_min)/x_gap
    y_max=np.max(y)
    y_min=np.min(y)
    y_gap=y_max-y_min
    if (y_gap==0):
        y_gap=1
    y=(y-y_min)/y_gap
    
    n=len(y)
    xy=np.dot(x,y)
    xx=np.dot(x,x)
    yy=np.dot(y,y)
    s_x=np.sum(x)
    s_y=np.sum(y)
    denominator = math.sqrt((n*xx-s_x*s_x)*(n*yy-s_y*s_y))
    if denominator == 0:
        return 0
    cc=(n*xy-s_x*s_y)/denominator
    return cc 

@numba.njit
def dcc(x,y):   # cross correlation without mean subtraction, suggested by Ambjornsson
    """Cosine between x and y WITHOUT centring them. cc_g.

    Because the mean is never removed, this responds to the absolute level and to the relative
    heights of the peaks, not only to where they sit. Two label sets at the same positions but
    with different brightness give different values. Returns 0 when either vector is all zero.
    """
    n=len(y)
    x_contig = np.ascontiguousarray(x)
    y_contig = np.ascontiguousarray(y)
    x_dot = np.dot(x_contig,x_contig)
    y_dot = np.dot(y_contig,y_contig)
    denominator = math.sqrt(x_dot) * math.sqrt(y_dot)
    if denominator == 0:
        return 0.0
    cc=np.dot(x_contig,y_contig)/denominator
    return cc

@numba.njit
def pcc(x,y):   # Pearson cross correlation
    """Alias for corrcoef. Not used by the distributed pipeline."""
    cc=corrcoef(x,y)
    return cc

#@numba.njit
def cross_correlation_using_fft(x, y):
    """Circular cross-correlation of two equal-length signals via FFT, zero lag at the centre.

    Only used by compute_shift, which the distributed pipeline does not call.
    """
    f1 = fft(x)
    f2 = fft(np.flipud(y))
    cc = np.real(ifft(f1 * f2))
    return fftshift(cc)
 
# shift &lt; 0 means that y starts 'shift' time steps before x # shift &gt; 0 means that y starts 'shift' time steps after x
#@numba.njit
def compute_shift(x, y):
    """Lag, in samples, that best aligns y onto x. Negative means y starts earlier.

    Not used by the distributed pipeline; the alignment code searches shifts explicitly instead.
    """
    assert len(x) == len(y)
    c = cross_correlation_using_fft(x, y)
    assert len(c) == len(x)
    zero_index = int(len(x) / 2) - 1
    shift = zero_index - np.argmax(c)
    return shift

@numba.njit
def Dcc(xx,y): # Ambjornsson suggestion
    """Slide the short vector y along the long vector xx and return the OFFSET of the best dcc.

    Returns the offset only, not the score. Not used by the distributed pipeline.
    """
    L=len(xx)
    n=len(y)
    max=0
    shift=0
    for i in range(L-n):
        x=xx[i:n+i]
        x_dot = np.dot(x,x)
        y_dot = np.dot(y,y)
        denominator = math.sqrt(x_dot) * math.sqrt(y_dot)
        if denominator == 0:
            D = 0.0
        else:
            D=np.dot(x,y)/denominator
        if (D>max):
            max=D
            shift=i
    return shift 
# Read tiff files in a folder by descending file size
def Read_TIFF_folder(path):
    """Full paths of the .tif files in `path`, largest file first.

    Size order, not name order, so the biggest molecules are processed first.
    """
    if (path[-1]!='/'):
        path += "/"
    #code exmaple: files = os.listdir(os.getcwd())
    folder=sorted(os.listdir(path))
    folder.sort(key=lambda f: os.stat(path+f).st_size, reverse=True)
    tif_files=[]
    for filename in folder:
        #if(filename.find('.tif') > 0):                    
        if(filename[-4:]=='.tif'):
            tif_files.append(path+filename)
    return tif_files

def _is_multipage_grayscale_mode(mode):
    """True for the PIL modes a 2-frame experimental TIFF reports, as opposed to an RGB render."""
    # Modes for experimental multipage TIFFs (16-bit or 8-bit)
    return mode in ('I;16B', 'I;16', 'I;16L', 'I', 'L', 'P')

def tif2RGimage(file):
    """Load a 2-frame TIFF as one RGB PIL image, plus its width and height.

    Frame 0 becomes red and frame 1 green; blue stays empty. Each channel is min-max normalised
    and scaled to 280 then clipped at 255, so the brightest few percent saturate and the rest of
    the range is visible. Green is thresholded first, exactly as in tif2RGprofile. For display
    only -- nothing is scored off this image.
    """
    image=Image.open(file)
    w,h=image.size
    if _is_multipage_grayscale_mode(image.mode):
        data=tiff.imread(file)
        r=data[0].astype(np.float64).flatten()
        r=r-r.min()
        r=normalize(r)
        r=r*280
        g=data[1].astype(np.float64).flatten()
        g=g-g.min()
        g=normalize(g)
        thres=g.mean()+(g.max()-g.mean())*0.1
        g[g<thres]=0
        g=g*280
        b=np.zeros(h*w)
        r[r>255]=255
        g[g>255]=255
        pixel=np.column_stack((r,g,b)).astype('uint8')
        image=Image.frombuffer('RGB', (w, h), pixel, 'raw', 'RGB', 0, 0)
    return image,w,h

# Make signal profiles from image
def tif2RGprofile(file):
    """Collapse a molecule or reference TIFF to (red profile, green profile, green peak indices).

    One value per image column, taken as the column maximum, so a molecule of width w gives two
    length-w arrays. Green is thresholded at mean + 0.1*(max - mean) and anything below is set to
    zero; the peaks returned are scipy find_peaks on that thresholded profile. Accepts both the
    2-frame experimental TIFFs and the RGB renders built from FASTA.
    """
    img=Image.open(file)
    if _is_multipage_grayscale_mode(img.mode): #Experimental Tiff files (8-bit or 16-bit multipage)
        data=tiff.imread(file).astype(np.float64)
        r=np.max(data[0], axis=0)
        g=np.max(data[1], axis=0)
        thres=g.mean()+(g.max()-g.mean())*0.1
        g[g<thres]=0
    elif (img.mode=='RGB'): #Image from FASTA file
        w,h=img.size
        color=img.split()
        r_data=np.array(color[0].getdata())
        r_data=r_data.reshape(h,w)
        r=np.max(r_data,axis=0)
        g_data=np.array(color[1].getdata())
        g_data=g_data.reshape(h,w)
        g=np.max(g_data,axis=0)
        thres=g.mean()+(g.max()-g.mean())*0.1
        g[g<thres]=0
    peaks,_= signal.find_peaks(g)
    return r,g, peaks

# Calculate Pearson Cross-correlation Value
@numba.njit
def crosscoef(ar1, ar2, t):
    """Pearson of ar1 and ar2 at lag t, with the number of overlapping samples.

    Positive t slides ar2 to the right. Not used by the distributed pipeline.
    """
    if (t>=0):
        x=ar1[t:min(len(ar1),len(ar2)+t)]
        y=ar2[0:min(len(ar2),len(ar1)-t)]
    elif(t<0):
        x=ar1[0:min(len(ar1),len(ar2)+t)]
        y=ar2[-t:min(len(ar2),len(ar1)-t)]
    cc=np.corrcoef(x,y)
    return cc[0,1],len(x)

#@numba.njit
def low_pass_filter(sg, order, cutoff):
    """Zero-phase Butterworth low-pass (filtfilt). Not used by the distributed pipeline."""
    b, a = signal.butter(order, cutoff) 
    lpf_sig = signal.filtfilt(b, a, sg) 
    return lpf_sig

def normalize(ar):
    """Min-max rescale to [0, 1]. A constant array comes back as zeros; an empty one is returned as is."""
    if (len(ar)==0):
        return ar
    min=ar.min()
    max=ar.max()
    if (max==min):
        max=max+1
    rt=(ar-min)/(max-min)
    return rt

def show_image(img):
    """Show one image full-bleed. Interactive helper, not used by the distributed pipeline."""
    fig=plt.figure(figsize=(10, 10))
    plt.imshow(img, cmap='gray')
    plt.axis('off')
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)    
    plt.show()

def show_image2(img1,img2):
    """Two images side by side. Interactive helper, not used by the distributed pipeline."""
    fig=plt.figure(figsize=(20, 10))
    plt.subplot(121)
    plt.axis('off')    
    plt.imshow(img1,cmap='gray')
    plt.subplot(122)
    plt.axis('off')
    plt.imshow(img2, cmap='gray')
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.show()
def show_image3(img1,img2,img3):
    """Three images in a row. Interactive helper, not used by the distributed pipeline."""
    fig=plt.figure(figsize=(20, 10))
    plt.subplot(131)
    plt.axis('off')    
    plt.imshow(img1,cmap='gray')
    plt.subplot(132)
    plt.axis('off')
    plt.imshow(img2, cmap='gray')
    plt.subplot(133)
    plt.axis('off')
    plt.imshow(img3, cmap='gray')
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.show()
def show_image4(img1,img2,img3,img4):
    """Four images in a 2x2. Interactive helper, not used by the distributed pipeline."""
    fig=plt.figure(figsize=(12,12))
    ax = fig.subplots(ncols=2, nrows=2)
    ax[0,0].axis('off')
    ax[0,0].imshow(img1, cmap='gray')
    ax[0,1].axis('off')
    ax[0,1].imshow(img2, cmap='gray')
    ax[1,0].axis('off')
    ax[1,0].imshow(img3, cmap='gray')
    ax[1,1].axis('off')
    ax[1,1].imshow(img4, cmap='gray')
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.show()    
# Draw reference and query image of aligned regions
def two_images(image1,w1,h1,image2,w2,h2, shift):
    """Stack two images vertically with a 10 px gap, offsetting the second by `shift` px.

    Not used by the distributed pipeline; the HTML view composes its own panels.
    """
    if (shift>=0):
        w=np.maximum(w1,w2+shift)
        x1=0
        x2=shift
    else:
        w=np.maximum(w1-shift, w2)
        x1=-shift
        x2=0
    two_image=Image.new("RGB",(w,h1+h2+10))
    two_image.paste(image1,(x1,0))
    two_image.paste(image2,(x2,h1+10))
    return two_image

def read_maps_file(map_file):
    """Split the FIRST line of a .maps file on tabs. Not used by the distributed pipeline."""
    f=open(map_file, encoding='utf-8')
    line=f.readline()
    frags=line.split("\t")
    return frags

def peaks2maps(tif_filename, DNA_size, peaks, map_filename):
    """Write one .maps line: name, DNA size, fragment count, then the inter-label fragment lengths.

    This is the only form maligner_dp reads. `peaks` and `DNA_size` must both already be in bp.
    n labels give n+1 fragments: start to the first label, each gap, and the last label to the end.
    """
    map_file=open(map_filename, "w", encoding='utf-8')
    frags=[]
    frags.append(tif_filename)
    frags.append(DNA_size)
    # Number of fragments = number of peaks + 1 (before first peak + between peaks + after last peak)
    num_fragments = len(peaks) + 1
    frags.append(num_fragments)
    
    # First fragment: from start (0) to first peak
    if len(peaks) > 0:
        frags.append(peaks[0])
    
    # Middle fragments: between consecutive peaks
    for i in range(len(peaks)-1):
        frags.append(peaks[i+1]-peaks[i])
    
    # Last fragment: from last peak to end (DNA_size)
    # Note: peaks are in pixel units, DNA_size is in bp units
    # peaks and DNA_size are both in bp, so the trailing fragment is DNA_size - peaks[-1].
    if len(peaks) > 0:
        # Calculate last fragment: DNA_size - last peak position
        last_fragment = DNA_size - peaks[-1]
        frags.append(last_fragment)
    
    str1=('\t'.join([str(x) for x in frags]))
    print (str1, file=map_file)
    map_file.close()

# select a ref map or make a new ref map
def use_tif_RGmap(map_tmp, quiet=False):
    """Return the reference map path, given either a .tif (used as is) or a .fasta.

    The .fasta branch is legacy: it shells out to a bare `seq_to_RGmap.py` with no arguments and
    is not the supported route. Build references with make_reference.py instead.
    """
    if not quiet:
        print(map_tmp)
    if map_tmp[-4:]=='.tif' and os.path.isfile(map_tmp):
        RGmap = map_tmp
    elif map_tmp[-6:]=='.fasta':  
        RGmap = map_tmp.replace(".fasta","_ref_RGmap.tif")
        if (os.path.exists(RGmap)==False):
            # making Reference map from fasta to tif
            print("Making a new reference map...\n")
            cmd = "seq_to_RGmap.py"
            os.system(cmd)
    else:
        print("Please input an existing fasta or tif file for a reference RGmap.")
        exit()

    if not quiet:
        print("Use:", RGmap)
    return RGmap