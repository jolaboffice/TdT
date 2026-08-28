#!/usr/bin/env python3
"""
hDOM_tif2html.py
Creates optimized HTML reports for hDOM analysis results and a summary index page.
Handles reverse matching direction (dir = -1).
"""

import os
import argparse
import io
import base64
import re
import shutil
import pandas as pd
import numpy as np
from PIL import Image, ImageEnhance
from DOM_lib_Q import tif2RGimage

# Constants
BPP = 200.0
BRIGHTNESS_ENHANCE = 1.5
REF_BRIGHTNESS_ENHANCE = 1.2

def image_to_base64(image):
    """Convert PIL Image to base64 string."""
    buf = io.BytesIO()
    image.save(buf, format='PNG')
    return base64.b64encode(buf.getvalue()).decode('utf-8')

def process_query_image(tif_path):
    """Load, enhance, and pad the query TIF image."""
    if not os.path.exists(tif_path):
        print(f"Warning: .tif file not found: {tif_path}")
        return None, None
    try:
        img, _, _ = tif2RGimage(tif_path)
        img = ImageEnhance.Brightness(img).enhance(BRIGHTNESS_ENHANCE)
        
        # Add 10px black padding
        w, h = img.size
        padded_w = w + 20
        new_img = Image.new('RGB', (padded_w, h), (0, 0, 0))
        new_img.paste(img, (10, 0))
        
        return new_img, image_to_base64(new_img)
    except Exception as e:
        print(f"Error processing query image: {e}")
        return None, None

def process_ref_match(match_info, ref_map_folder, query_final_width):
    """Crop, enhance, and flip reference matching region based on direction."""
    ref_name = match_info['ref_name']
    aligned_bp = match_info['aligned_bp']
    scale = match_info['scale']
    direction = match_info.get('dir', 1)
    
    ref_path = os.path.join(ref_map_folder, f"{ref_name}.tif")
    if not os.path.exists(ref_path):
        return None

    try:
        ref_img, ref_w, ref_h = tif2RGimage(ref_path)
        shift_px = int(aligned_bp / BPP)
        
        # Calculate region width
        query_w = max(query_final_width - 20, 1)
        ref_region_w = int(query_w * scale)
        
        # Crop with 10px buffer
        start_x = max(0, shift_px - 10)
        end_x = min(ref_w, shift_px + ref_region_w + 10)
        
        if start_x >= end_x:
            return None
            
        crop = ref_img.crop((start_x, 0, end_x, ref_h))
        
        # Handle reverse direction
        if direction == -1:
            crop = crop.transpose(Image.FLIP_LEFT_RIGHT)
            
        crop = crop.resize((query_final_width, ref_h))
        
        # Channel adjustments
        arr = np.array(crop)
        
        # Normalize Red
        red = arr[:, :, 0].astype(np.float32)
        r_min, r_max = np.min(red), np.max(red)
        if r_max > r_min:
            red = ((red - r_min) / (r_max - r_min)) * 255.0
        arr[:, :, 0] = np.clip(red, 0, 255).astype(np.uint8)
        
        # Increase Green
        arr[:, :, 1] = np.clip((arr[:, :, 1].astype(np.float32) * 2.0), 0, 255).astype(np.uint8)
        
        enhanced = Image.fromarray(arr)
        enhanced = ImageEnhance.Brightness(enhanced).enhance(REF_BRIGHTNESS_ENHANCE)
        
        return image_to_base64(enhanced)
    except Exception as e:
        print(f"Error processing ref {ref_name}: {e}")
        return None

def get_table_html(df, folder_name="", file_name=""):
    """Format and generate HTML table from matching results with hyperlinks on rank."""
    if df is None or df.empty:
        return "<p>No matching results available.</p>"
        
    df_fmt = df.copy()
    cols = df_fmt.columns
    
    # Formatting
    if 'aligned_bp' in cols:
        df_fmt['aligned_bp'] = df_fmt['aligned_bp'].apply(lambda x: f"{int(x):,}" if pd.notna(x) else "")
    if 'shift_bp' in cols:
        df_fmt['shift_bp'] = df_fmt['shift_bp'].apply(lambda x: f"{int(x):,}" if pd.notna(x) else "")
    
    for c in ['scale', 'cc_r', 'cc_g', 'log(P)']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(lambda x: f"{x:.2f}" if pd.notna(x) and isinstance(x, (int, float)) else x)
            
    for c in ['cc_rg', 'cc_rg2']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(lambda x: f"{x:.1e}" if pd.notna(x) and isinstance(x, (int, float)) else x)
            
    # Hyperlink for rank. 'sel' stays the table row position (1,2,...) so view.php
    # still selects the right row, but the visible label uses the per-method rank:
    # PCC rows -> 1,2,...  maligner rows -> m1,m2,...
    if 'rank' in cols:
        methods = df['method'].tolist() if 'method' in df.columns else [''] * len(df_fmt)
        ranks = df['rank'].tolist() if 'rank' in df.columns else list(range(1, len(df_fmt) + 1))
        new_ranks = []
        for i in range(len(df_fmt)):
            display_rank = i + 1  # row position for view.php sel
            method = str(methods[i]).lower() if i < len(methods) else ''
            try:
                rk = int(ranks[i])
            except (TypeError, ValueError):
                rk = display_rank
            label = f"m{rk}" if 'mal' in method else f"{rk}"
            new_ranks.append(f'<a href="../view.php?folder={folder_name}&file={file_name}&sel={display_rank}">{label}</a>')
        df_fmt['rank'] = new_ranks

    for c in ['misc', 'dir']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(lambda x: f"{int(x)}" if pd.notna(x) else "")

    # Display-only header relabels: log(P) -> "score", ref_name -> "chromosome"
    df_fmt = df_fmt.rename(columns={'log(P)': 'score', 'ref_name': 'chromosome'})

    html = df_fmt.to_html(index=False, classes='matching-table', border=0, escape=False)
    
    # Add conditional row classes (PCC/MAL coloring and Rank highlighting)
    if 'method' in cols:
        parts = html.split('<tbody>')
        if len(parts) == 2:
            header, body = parts[0] + '<tbody>', parts[1]
            rows = re.findall(r'<tr>.*?</tr>', body, re.DOTALL)
            new_rows = []
            for i, row in enumerate(rows):
                rank_raw = df.iloc[i].get('rank', 0)
                m_val = str(df.iloc[i]['method']).lower() if pd.notna(df.iloc[i]['method']) else ''
                
                cls = []
                if 'pcc' in m_val: cls.append('method-pcc')
                elif 'mal' in m_val: cls.append('method-mal')
                
                try:
                    # Highlight Top 2 (Rank 1 and 2 in CSV)
                    if int(rank_raw) in (1, 2):
                        cls.append('rank-top')
                except (TypeError, ValueError):
                    pass
                
                if cls:
                    row = row.replace('<tr>', f'<tr class="{" ".join(cls)}">', 1)
                new_rows.append(row)
            
            # Reconstruct with closing tags
            footer_match = re.search(r'(</tbody>.*?</table>)', body, re.DOTALL)
            footer = footer_match.group(1) if footer_match else '</tbody></table>'
            html = header + "".join(new_rows) + footer
            
    return html

def get_html_style():
    """Return the CSS styles for the report."""
    return """
    body { font-family: 'Segoe UI', sans-serif; line-height: 1.6; color: #333; max-width: 1200px; margin: 0 auto; padding: 40px 20px; background-color: #f9f9f9; }
    .container { background-color: #fff; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #2980b9; margin-top: 30px; }
    .info-box { background-color: #f8f9fa; border-left: 5px solid #3498db; padding: 15px; margin-bottom: 20px; white-space: pre-wrap; font-family: monospace; font-size: 0.95em; overflow-x: auto; }
    .image-container { margin: 30px 0; text-align: center; background-color: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd; }
    .comparison-box { margin: 30px 0; padding: 20px; background-color: #000; border-radius: 5px; }
    .comparison-images { display: flex; flex-direction: column; gap: 0; }
    .comparison-item { display: flex; align-items: center; width: 100%; }
    .image-label { min-width: 60px; padding: 0 15px; font-weight: bold; color: #fff; text-align: right; }
    .comparison-item .image-container { margin: 0; padding: 0; border: none; flex: 1; background-color: #000; }
    .comparison-item img { display: block; margin: 0 auto; max-width: 100%; height: auto; }
    .table-container { overflow-x: auto; margin-top: 15px; }
    .matching-table { width: 100%; border-collapse: collapse; font-size: 0.9em; }
    .matching-table th { background-color: #3498db; color: #fff; text-align: right; padding: 12px; position: sticky; top: 0; }
    .matching-table th:first-child { text-align: left; }
    .matching-table td { padding: 10px 12px; border-bottom: 1px solid #eee; text-align: right; }
    .matching-table td:first-child { text-align: left; }
    .matching-table tr:hover { background-color: #f1f7ff; }
    .method-pcc { background-color: #e8f4f8; }
    .method-pcc:nth-child(even) { background-color: #d4e9f0; }
    .method-mal { background-color: #fff4e6; }
    .method-mal:nth-child(even) { background-color: #ffe8cc; }
    .rank-top { font-weight: bold; }
    .summary-table { width: 100%; border-collapse: collapse; margin-top: 20px; }
    .summary-table th { background-color: #2c3e50; color: #fff; padding: 12px; text-align: left; cursor: pointer; user-select: none; }
    .summary-table td { padding: 10px 12px; border-bottom: 1px solid #ddd; }
    .summary-table tr:hover { background-color: #f5f5f5; }
    .text-right { text-align: right !important; }
    .text-center { text-align: center !important; }
    footer { margin-top: 50px; font-size: 0.8em; color: #888; text-align: center; }
    """

def tif2html(tif_path):
    """Generates an HTML report and returns summary information."""
    ext = os.path.splitext(tif_path)[1].lower()
    if ext not in (".tif", ".tiff"):
        print(f"Error: Input must be a .tif or .tiff file. Got: {tif_path}")
        return None
    
    base_path = os.path.splitext(tif_path)[0]
    info_path = f"{base_path}_info.txt"
    csv_path = f"{base_path}_matching_results.csv"
    html_path = f"{base_path}.html"

    # 1. Load Info
    info_content = f"Information file not found: {info_path}"
    if os.path.exists(info_path):
        with open(info_path, 'r', encoding='utf-8') as f:
            info_content = f.read()

    # 2. Process Images
    query_img, query_b64 = process_query_image(tif_path)
    if not query_img:
        return None
    
    query_final_w = query_img.size[0]
    molecule_length_bp = int((query_final_w - 20) * BPP)
    
    # 3. Read Results & Generate Comparisons
    df = pd.read_csv(csv_path) if os.path.exists(csv_path) else None
    comp_html = ""
    summary_info = {
        'tif_name': os.path.basename(tif_path),
        'html_link': os.path.basename(html_path),
        'ref_map': 'N/A',
        'aligned_bp': 'N/A',
        'length_bp': molecule_length_bp,
        'direction': '+',
        'log_p': 'N/A'
    }
    
    if df is not None and not df.empty:
        top_hit = df.iloc[0]
        summary_info['ref_map'] = top_hit.get('ref_name', 'N/A')
        summary_info['log_p'] = top_hit.get('log(P)', 'N/A')
        aligned = top_hit.get('aligned_bp')
        if pd.isna(aligned):
            aligned = top_hit.get('shift_bp', 0)
        summary_info['aligned_bp'] = int(aligned)
        
        # Handle direction in summary
        direction = top_hit.get('dir', 1)
        summary_info['direction'] = '+' if direction == 1 else '-'

        # Find ref_map folder
        tif_dir = os.path.dirname(os.path.abspath(tif_path))
        ref_map_folder = os.path.join(os.path.dirname(tif_dir) or os.getcwd(), "ref_map")
        if not os.path.isdir(ref_map_folder): ref_map_folder = "ref_map"
        
        # Prepare Top 2 hits
        comp_rows = []
        for i in range(min(2, len(df))):
            row = df.iloc[i]
            ref_name = row.get('ref_name')
            if pd.isna(ref_name) or not ref_name:
                continue
            aligned_bp = row.get('aligned_bp')
            if pd.isna(aligned_bp):
                aligned_bp = row.get('shift_bp', 0)
            scale = row.get('scale', 1.0)
            if pd.isna(scale):
                scale = 1.0
            dir_val = row.get('dir', 1)
            comp_rows.append({'ref_name': ref_name, 'aligned_bp': aligned_bp, 'scale': scale, 'dir': dir_val})
            
        ref_htmls = ""
        for i, row in enumerate(comp_rows):
            ref_b64 = process_ref_match(row, ref_map_folder, query_final_w)
            if ref_b64:
                ref_htmls += f'<div class="comparison-item"><div class="image-label">rank{i+1}</div>' \
                             f'<div class="image-container"><img src="data:image/png;base64,{ref_b64}"></div></div>'
        
        comp_html = f'<div class="comparison-box"><div class="comparison-images">' \
                    f'<div class="comparison-item"><div class="image-label">{molecule_length_bp/1000:.0f} kb</div>' \
                    f'<div class="image-container"><img src="data:image/png;base64,{query_b64}"></div></div>' \
                    f'{ref_htmls}</div></div>'
    else:
        comp_html = f'<div class="image-container"><img src="data:image/png;base64,{query_b64}"></div>'

    # 4. Generate Final HTML
    tif_abs = os.path.abspath(tif_path)
    folder_name = os.path.basename(os.path.dirname(tif_abs))
    file_name = os.path.basename(tif_abs)
    table_html = get_table_html(df, folder_name=folder_name, file_name=file_name)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>hDOM Analysis - {os.path.basename(tif_path)}</title>
    <style>{get_html_style()}</style>
</head>
<body>
    <div class="container">
        <h1>hDOM Analysis Report</h1>
        <p><strong>File:</strong> {os.path.basename(tif_path)}</p>
        {comp_html}
        <div class="table-container">{table_html}</div>
        <h2>Information</h2>
        <div class="info-box">{info_content}</div>
    </div>
    <footer>Generated by hDOM_tif2html.py</footer>
</body>
</html>"""

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created: {html_path}")
    return summary_info

def generate_index_html(summaries, output_path):
    """Creates a summary index HTML page with sortable columns."""
    rows_html = ""
    for s in summaries:
        aligned_bp_fmt = f"{s['aligned_bp']:,}" if isinstance(s['aligned_bp'], (int, float)) else s['aligned_bp']
        length_bp_fmt = f"{s['length_bp']:,}" if isinstance(s['length_bp'], (int, float)) else s['length_bp']
        log_p_fmt = f"{s['log_p']:.2f}" if isinstance(s['log_p'], (int, float)) else s['log_p']
        rows_html += f"""
        <tr>
            <td data-value="{s['tif_name'].lower()}"><a href="{s['html_link']}">{s['tif_name']}</a></td>
            <td data-value="{s['ref_map'].lower()}">{s['ref_map']}</td>
            <td class="text-right" data-value="{s['aligned_bp']}">{aligned_bp_fmt}</td>
            <td class="text-right" data-value="{s['length_bp']}">{length_bp_fmt}</td>
            <td class="text-right" data-value="{s['log_p']}">{log_p_fmt}</td>
            <td class="text-center" data-value="{s['direction']}">{s['direction']}</td>
        </tr>"""

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>hDOM Analysis Summary Index</title>
    <style>{get_html_style()}</style>
    <script>
    function sortTable(n) {{
        var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
        table = document.getElementById("summaryTable");
        switching = true;
        dir = "asc";
        while (switching) {{
            switching = false;
            rows = table.rows;
            for (i = 1; i < (rows.length - 1); i++) {{
                shouldSwitch = false;
                x = rows[i].getElementsByTagName("TD")[n];
                y = rows[i + 1].getElementsByTagName("TD")[n];
                
                var xVal = x.getAttribute("data-value");
                var yVal = y.getAttribute("data-value");
                
                // Multi-level sort for Aligned Position (column 2): 
                // Primary: Reference Map (column 1), Secondary: Aligned Position (column 2)
                if (n === 2) {{
                    var xRef = rows[i].getElementsByTagName("TD")[1].getAttribute("data-value").toLowerCase();
                    var yRef = rows[i + 1].getElementsByTagName("TD")[1].getAttribute("data-value").toLowerCase();
                    if (xRef !== yRef) {{
                        if (dir == "asc") {{
                            if (xRef > yRef) {{ shouldSwitch = true; break; }}
                        }} else {{
                            if (xRef < yRef) {{ shouldSwitch = true; break; }}
                        }}
                        continue; 
                    }}
                }}
                
                // Try numeric comparison if both are numbers
                var xNum = parseFloat(xVal);
                var yNum = parseFloat(yVal);
                
                if (!isNaN(xNum) && !isNaN(yNum)) {{
                    if (dir == "asc") {{
                        if (xNum > yNum) {{ shouldSwitch = true; break; }}
                    }} else if (dir == "desc") {{
                        if (xNum < yNum) {{ shouldSwitch = true; break; }}
                    }}
                }} else {{
                    if (dir == "asc") {{
                        if (xVal.toLowerCase() > yVal.toLowerCase()) {{ shouldSwitch = true; break; }}
                    }} else if (dir == "desc") {{
                        if (xVal.toLowerCase() < yVal.toLowerCase()) {{ shouldSwitch = true; break; }}
                    }}
                }}
            }}
            if (shouldSwitch) {{
                rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                switching = true;
                switchcount ++;
            }} else {{
                if (switchcount == 0 && dir == "asc") {{
                    dir = "desc";
                    switching = true;
                }}
            }}
        }}
    }}
    </script>
</head>
<body>
    <div class="container">
        <h1>hDOM Analysis Summary Index</h1>
        <table class="summary-table" id="summaryTable">
            <thead>
                <tr>
                    <th onclick="sortTable(0)">TIF File</th>
                    <th onclick="sortTable(1)">Reference Map</th>
                    <th class="text-right" onclick="sortTable(2)">Aligned Position (bp)</th>
                    <th class="text-right" onclick="sortTable(3)">Molecule Length (bp)</th>
                    <th class="text-right" onclick="sortTable(4)">score</th>
                    <th class="text-center" onclick="sortTable(5)">Dir</th>
                </tr>
            </thead>
            <tbody>
                {rows_html}
            </tbody>
        </table>
    </div>
    <footer>Generated by hDOM_tif2html.py</footer>
</body>
</html>"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created Summary Index: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Convert hDOM .tif results to HTML reports.")
    parser.add_argument("tif_files", nargs='+', help="One or more .tif files or a directory containing .tif files")
    args = parser.parse_args()
    
    all_tif_files = []
    for path in args.tif_files:
        if os.path.isdir(path):
            files = [os.path.join(path, f) for f in os.listdir(path) if f.lower().endswith(('.tif', '.tiff'))]
            all_tif_files.extend(sorted(files))
        else:
            all_tif_files.append(path)
    
    if not all_tif_files:
        print("No .tif files found to process.")
        return

    summaries = []
    for f in all_tif_files:
        summary = tif2html(f)
        if summary:
            summaries.append(summary)
    
    if summaries:
        # Generate index.html in the same directory as the first processed file or current dir
        first_file_dir = os.path.dirname(os.path.abspath(all_tif_files[0]))
        index_path = os.path.join(first_file_dir, "index.html")
        generate_index_html(summaries, index_path)

if __name__ == "__main__":
    main()
