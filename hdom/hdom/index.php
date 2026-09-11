<?php
register_shutdown_function(function(){ if(($_SERVER['REQUEST_METHOD']??'')==='POST') return; $h=(($_GET['file']??'')!=='')?'index.php':'../'; echo '<a href="'.$h.'" style="position:fixed;top:6px;left:6px;z-index:2147483647;background:#34495e;color:#fff;padding:6px 14px;border-radius:5px;text-decoration:none;font:bold 13px sans-serif;box-shadow:0 1px 4px rgba(0,0,0,.4)">&#8593; Up</a>'; });
// 4-method local-align comparison index. Serves <base>_view_ytb.html, or <base>_view.html:
// the analysis server writes the first, hDOM_pipeline.py writes the second, and one viewer has to
// read both or a released run lands on a 404.
// and re-renders the matching table (local / mDOM / mal / pcc + H) from the CSV.
$base_path = __DIR__;
header('Cache-Control: no-store, no-cache, must-revalidate');

// ---- mark-root: all *_marks.json live at the DOM_local root (the ancestor dir named 'local'),
// NOT dirname($base_path). Nested folders (zave/fresh_reclass/A_to_C, zave/C_total/*, ...) are
// several levels below 'local'; without this their marks would try to write into their own parent
// (non-existent, apache-unwritable) json and silently vanish, and reclassify (which reads local/*)
// would never see them. Walk up to 'local' so every folder shares the same canonical mark files.
$__p = $base_path;
while ($__p && $__p !== '/' && $__p !== '.' && basename($__p) !== 'local' && basename($__p) !== 'ML_label_local_ccr') { $__p = dirname($__p); }
$MARK_ROOT = (basename($__p) === 'local' || basename($__p) === 'ML_label_local_ccr') ? $__p : dirname($base_path);

// ---- Cplus marking (C group only): server-side flag, survives navigation ----
// Marks live at local/cplus_marks.json (root, survives reclassify wipe) and must be
// pre-created world-writable (666) by reclassify so apache can overwrite them.
// C-marking enabled for the C group, C_total subfolders, and ANY fresh_reclass transition bucket
// (buckets are per-transition like C_to_B / A_to_C, so match by parent dir rather than listing each).
$IS_C = in_array(basename($base_path), array('C', 'C_longmolecules', 'C_longer_than_200kb', 'candidate_recovered', 'linearization', 'others', 'Long_Molecules', 'C_freshscan_income'), true)
     || basename(dirname($base_path)) === 'fresh_reclass';
$CPLUS_JSON = $MARK_ROOT . '/cplus_marks.json';
function cplus_load($p){ if(!is_file($p)) return array(); $j=json_decode((string)file_get_contents($p),true); return is_array($j)?$j:array(); }
function cplus_save($p,$j){ $fp=fopen($p,'c+'); if(!$fp) return; flock($fp,LOCK_EX); ftruncate($fp,0); rewind($fp);
    fwrite($fp,json_encode($j)); fflush($fp); flock($fp,LOCK_UN); fclose($fp); }
// read one group's _summary.csv -> array(molecule => full row line array), plus header
function summ_read($csv){ $rows=array(); $head=null; if(($h=@fopen($csv,'r'))){ $head=fgetcsv($h);
    while(($r=fgetcsv($h))!==false){ if($head && !empty($r[0])) $rows[$r[0]]=$r; } fclose($h); } return array($head,$rows); }
function summ_write($csv,$head,$rows){ if(($h=fopen($csv,'w'))){ if($head)fputcsv($h,$head);
    foreach($rows as $r) fputcsv($h,$r); fclose($h); } }

if ($_SERVER['REQUEST_METHOD']==='POST' && $IS_C && ($_POST['action']??'')==='toggle_cplus') {
    header('Content-Type: application/json');
    $name=(string)($_POST['name']??'');
    if($name===''||strpbrk($name,"/\\")!==false||strpos($name,'..')!==false){ echo json_encode(['ok'=>false]); exit; }
    $j=cplus_load($CPLUS_JSON);
    if(!empty($j[$name])) unset($j[$name]); else $j[$name]=true;
    cplus_save($CPLUS_JSON,$j);
    echo json_encode(['ok'=>true,'marked'=>!empty($j[$name]),'count'=>count($j)]); exit;
}
// ---- 't' trash/discard marking (all groups): root trash_marks.json, survives navigation & reclassify ----
$TRASH_JSON = $MARK_ROOT . '/trash_marks.json';
if ($_SERVER['REQUEST_METHOD']==='POST' && ($_POST['action']??'')==='toggle_trash') {
    header('Content-Type: application/json');
    $name=(string)($_POST['name']??'');
    if($name===''||strpbrk($name,"/\\")!==false||strpos($name,'..')!==false){ echo json_encode(['ok'=>false]); exit; }
    $j=cplus_load($TRASH_JSON);
    if(!empty($j[$name])) unset($j[$name]); else $j[$name]=true;
    cplus_save($TRASH_JSON,$j);
    echo json_encode(['ok'=>true,'trashed'=>!empty($j[$name]),'count'=>count($j)]); exit;
}
// ---- 'm' mutation/suspicious marking (all groups): root m_marks.json ----
$M_JSON = $MARK_ROOT . '/m_marks.json';
if ($_SERVER['REQUEST_METHOD']==='POST' && ($_POST['action']??'')==='toggle_m') {
    header('Content-Type: application/json');
    $name=(string)($_POST['name']??'');
    if($name===''||strpbrk($name,"/\\")!==false||strpos($name,'..')!==false){ echo json_encode(['ok'=>false]); exit; }
    $j=cplus_load($M_JSON);
    if(!empty($j[$name])) unset($j[$name]); else $j[$name]=true;
    cplus_save($M_JSON,$j);
    echo json_encode(['ok'=>true,'marked'=>!empty($j[$name]),'count'=>count($j)]); exit;
}
// ---- 'a' → A marking (all groups): root a_marks.json (move to A group manually) ----
$A_JSON = $MARK_ROOT . '/a_marks.json';
if ($_SERVER['REQUEST_METHOD']==='POST' && ($_POST['action']??'')==='toggle_a') {
    header('Content-Type: application/json');
    $name=(string)($_POST['name']??'');
    if($name===''||strpbrk($name,"/\\")!==false||strpos($name,'..')!==false){ echo json_encode(['ok'=>false]); exit; }
    $j=cplus_load($A_JSON);
    if(!empty($j[$name])) unset($j[$name]); else $j[$name]=true;
    cplus_save($A_JSON,$j);
    echo json_encode(['ok'=>true,'marked'=>!empty($j[$name]),'count'=>count($j)]); exit;
}
// ---- 'k' → kink marking (all groups): root kink_marks.json (curvy/bent molecules) ----
$K_JSON = $MARK_ROOT . '/kink_marks.json';
if ($_SERVER['REQUEST_METHOD']==='POST' && ($_POST['action']??'')==='toggle_k') {
    header('Content-Type: application/json');
    $name=(string)($_POST['name']??'');
    if($name===''||strpbrk($name,"/\\")!==false||strpos($name,'..')!==false){ echo json_encode(['ok'=>false]); exit; }
    $j=cplus_load($K_JSON);
    if(!empty($j[$name])) unset($j[$name]); else $j[$name]=true;
    cplus_save($K_JSON,$j);
    echo json_encode(['ok'=>true,'marked'=>!empty($j[$name]),'count'=>count($j)]); exit;
}
// ---- 'r' → AR marking (all groups): root r_marks.json (acrocentric p-arm / rDNA region) ----
$R_JSON = $MARK_ROOT . '/r_marks.json';
if ($_SERVER['REQUEST_METHOD']==='POST' && ($_POST['action']??'')==='toggle_r') {
    header('Content-Type: application/json');
    $name=(string)($_POST['name']??'');
    if($name===''||strpbrk($name,"/\\")!==false||strpos($name,'..')!==false){ echo json_encode(['ok'=>false]); exit; }
    $j=cplus_load($R_JSON);
    if(!empty($j[$name])) unset($j[$name]); else $j[$name]=true;
    cplus_save($R_JSON,$j);
    echo json_encode(['ok'=>true,'marked'=>!empty($j[$name]),'count'=>count($j)]); exit;
}
// NOTE: the actual C->Cplus move is done by move_cplus.py (run as the file owner), because
// apache cannot create/delete files in these triplecells-owned dirs. The web UI only marks.

function hdom_table_from_csv($csv) {
    $fh = @fopen($csv, 'r'); if (!$fh) return null;
    $head = fgetcsv($fh); if (!$head) { fclose($fh); return null; }
    $rows = array();
    while (($r = fgetcsv($fh)) !== false) {
        if (count($r) < count($head)) $r = array_pad($r, count($head), '');
        $rows[] = array_combine($head, array_slice($r, 0, count($head)));
    }
    fclose($fh); if (!$rows) return null;

    $f2 = function($v){ return is_numeric($v) ? sprintf('%.2f', (float)$v) : htmlspecialchars((string)$v); };
    $fe = function($v){ return is_numeric($v) ? sprintf('%.1e', (float)$v) : htmlspecialchars((string)$v); };

    // group (keep original index for data-sel / click-override)
    $loc=array(); $mal=array();
    foreach ($rows as $i=>$r) {
        $mth = strtolower((string)($r['method'] ?? ''));
        if     (strpos($mth,'local')!==false) $loc[]=array($i,$r);
        elseif (strpos($mth,'mal')!==false)   $mal[]=array($i,$r);
    }
    $f2c = function($r,$k) use($f2){ return (isset($r[$k]) && $r[$k]!=='') ? $f2($r[$k]) : ''; };
    $fec = function($r,$k) use($fe){ return (isset($r[$k]) && $r[$k]!=='') ? $fe($r[$k]) : ''; };
    $bpf = function($r){ return (isset($r['aligned_bp']) && is_numeric($r['aligned_bp'])) ? number_format((int)$r['aligned_bp']) : htmlspecialchars((string)($r['aligned_bp']??'')); };
    $drf = function($r){ return (isset($r['dir']) && is_numeric($r['dir'])) ? (string)(int)$r['dir'] : htmlspecialchars((string)($r['dir']??'')); };
    $txf = function($r,$k){ return (isset($r[$k]) && $r[$k]!=='') ? htmlspecialchars((string)$r[$k]) : ''; };

    // Row colours are unchanged: hDOM=method-local(#e8f4f8), mal=method-mal(#fff4e6). Only the mal header is a muted tan.
    $css = '<style>'
      .'.mm .hdr-hdom td{background:#3498db;color:#fff;font-weight:bold}'
      .'.mm .hdr-mal td{background:#a67c52;color:#fff;font-weight:bold}'
      .'.mm .osep td{background:#34495e;color:#fff;font-weight:bold;text-align:left}'
      .'.mm .hdr-hdom,.mm .hdr-mal,.mm .osep{cursor:default}'
      .'.mm td:nth-child(n+6),.mm th:nth-child(n+6){width:9.5%;min-width:56px;text-align:center;box-sizing:border-box}'  // scale..h_score and the mal columns share one width (7 cols x 9.5% = 66.5%)
      .'.mm thead th,.mm .hdr-hdom td,.mm .hdr-mal td{padding-top:2px !important;padding-bottom:2px !important;line-height:1.25 !important}'  // header row: same font size, thinner vertically only
      // Data rows (hDOM/mal) only, at 0.9x height: padding 10+10 + line 17.3 = 37.3px -> 8.1+8.1+17.3 = 33.5px.
      // The header row (thead th), the group headers (hdr-*) and the Other ranks separator (osep) are
      // excluded, so the rule above still applies to them.
      .'.mm tbody tr:not(.hdr-hdom):not(.hdr-mal):not(.osep) td{padding-top:8.1px !important;padding-bottom:8.1px !important;line-height:1.2 !important}'
      .'</style>';

    // hDOM, 12 columns (cc_rg and misc dropped). mal leaves the cc_r slot (7th) empty and lines
    // matched..m_score up from there.
    $HDR_H = array('rank','chromosome','method','aligned_bp','dir','scale','cc_r','cc_g','cc_rg2','cc_rL','F1_g','h_score');
    $HDR_M = array('rank','chromosome','method','aligned_bp','dir','scale','matched','FP','FN','sizing','F1','m_score');
    $hrow = function($cols,$cls){ $t='    <tr class="'.$cls.'">'; foreach($cols as $c) $t.='<td>'.htmlspecialchars($c).'</td>'; return $t."</tr>\n"; };
    $rLoc = function($pair) use($f2c,$fec,$bpf,$drf){
        list($i,$r)=$pair; $rk=(string)($r['rank']??'');
        $cl='method-local'.((is_numeric($rk)&&((int)$rk===1||(int)$rk===2))?' rank-top':'');
        $tds='<td>rank'.htmlspecialchars(is_numeric($rk)?(string)(int)$rk:$rk).'</td>'
            .'<td>'.htmlspecialchars((string)($r['ref_name']??'')).'</td><td>hDOM</td>'
            .'<td>'.$bpf($r).'</td><td>'.$drf($r).'</td><td>'.$f2c($r,'scale').'</td>'
            .'<td>'.$f2c($r,'cc_r').'</td><td>'.$f2c($r,'cc_g').'</td><td>'.$fec($r,'cc_rg2').'</td>'
            .'<td>'.$f2c($r,'local_cc_r').'</td><td>'.$f2c($r,'F1').'</td><td>'.$f2c($r,'log(P)').'</td>';
        return '    <tr class="'.$cl.'" data-rtype="local" data-sel="'.($i+1).'">'.$tds."</tr>\n";
    };
    $rMal = function($pair) use($f2c,$txf,$bpf,$drf){
        list($i,$r)=$pair; $rk=(string)($r['rank']??'');
        $cl='method-mal'.((is_numeric($rk)&&((int)$rk===1||(int)$rk===2))?' rank-top':'');
        $mM=(int)($r['mal_matched']??0); $mFP=(int)($r['mal_FP']??0); $mFN=(int)($r['mal_FN']??0);
        $mden=2*$mM+$mFP+$mFN; $mf1=($mden>0)?number_format(2*$mM/$mden,2):'';   // F1 = 2M/(2M+FP+FN)
        $tds='<td>m'.htmlspecialchars(is_numeric($rk)?(string)(int)$rk:$rk).'</td>'
            .'<td>'.htmlspecialchars((string)($r['ref_name']??'')).'</td><td>mal</td>'
            .'<td>'.$bpf($r).'</td><td>'.$drf($r).'</td><td>'.$f2c($r,'scale').'</td>'
            .'<td>'.$txf($r,'mal_matched').'</td><td>'.$txf($r,'mal_FP').'</td><td>'.$txf($r,'mal_FN').'</td>'
            .'<td>'.$f2c($r,'mal_sizing').'</td><td>'.$mf1.'</td><td>'.$f2c($r,'m_score').'</td>';
        return '    <tr class="'.$cl.'" data-rtype="mal" data-sel="'.($i+1).'">'.$tds."</tr>\n";
    };

    // One table, 12 columns lining up on their own. thead is the hDOM header.
    $t = '<table border="0" class="dataframe matching-table mm">'."\n  <thead>\n    <tr>\n";
    foreach($HDR_H as $c) $t.='      <th>'.htmlspecialchars($c)."</th>\n"; $t.="    </tr>\n  </thead>\n  <tbody>\n";
    foreach (array_slice($loc,0,3) as $p) $t.=$rLoc($p);          // top: hDOM top-3
    $t.=$hrow($HDR_M,'hdr-mal');
    foreach (array_slice($mal,0,3) as $p) $t.=$rMal($p);          //      mal top-3
    $lr=array_slice($loc,3); $mr=array_slice($mal,3);
    if ($lr || $mr){                                              // below: Other ranks
        $t.='    <tr class="osep"><td colspan="'.count($HDR_H).'">Other ranks</td></tr>'."\n";
        if ($lr){ $t.=$hrow($HDR_H,'hdr-hdom'); foreach($lr as $p) $t.=$rLoc($p); }
        if ($mr){ $t.=$hrow($HDR_M,'hdr-mal'); foreach($mr as $p) $t.=$rMal($p); }
    }
    $t.="  </tbody>\n</table>";
    return $css . $t;
}

// ---- molecule view: index.php?file=<base> ----
if (isset($_GET['file']) && $_GET['file'] !== '') {
    $f = preg_replace('/\.(tiff?|html)$/i','', (string)$_GET['file']);
    $f = preg_replace('/_view(_ytb)?$/','',$f);
    if ($f==='' || strpbrk($f,"/\\")!==false || strpos($f,'..')!==false) { http_response_code(400); echo 'bad'; exit; }
    $vp = "$base_path/{$f}_view_ytb.html";
    if (!is_file($vp)) $vp = "$base_path/{$f}_view.html";
    if (!is_file($vp)) { http_response_code(404); echo 'view not found'; exit; }
    $m = file_get_contents($vp);
    $m = preg_replace('#<script data-vredirect>.*?</script>#s', '', $m);  // strip direct-access redirect when served here
    $m = preg_replace('#<script>\s*\(function\(\)\{(?:(?!</script>).)*?save_trash\.php.*?</script>#s', '', $m);  // strip baked (non-persisting) trash script; replaced below
    // fixed-width label column so data & all ref_map panels share left/right edges.
    // NEVER wrap (keep "mDOM1" intact); only the long H_local/H_mDOM labels get a smaller font.
    $cssfix = '<style>.comparison-item[data-rtype="mDOM"],.comparison-item[data-rtype="hDOM"],.comparison-item[data-rtype="pcc"]{display:none!important}.matching-table tr[data-rtype="mDOM"],.matching-table tr[data-rtype="hDOM"],.matching-table tr[data-rtype="pcc"]{display:none!important}.comparison-item .image-label{width:90px!important;min-width:90px!important;max-width:90px!important;'
            . 'box-sizing:border-box!important;white-space:nowrap!important;overflow:visible!important;}'
            . '.matching-table .method-local{background-color:#e8f4f8}.matching-table .method-local:nth-child(even){background-color:#d4e9f0}</style>';
    $m = preg_match('/<\/head>/i', $m) ? preg_replace('/<\/head>/i', $cssfix.'</head>', $m, 1) : $cssfix.$m;
    // shrink only the long homolog labels so they fit on one line (no wrap, no push)
    $m = preg_replace('/<div class="image-label">(H_[^<]*)<\/div>/',
                      '<div class="image-label" style="font-size:9px;padding:0 4px;letter-spacing:-.2px">$1</div>', $m);
    // ---- match reference layout: local1/2 -> rank1/2 panel labels; positional panel separators ----
    $m = str_replace(array('<div class="image-label">local1</div>','<div class="image-label">local2</div>'),
                     array('<div class="image-label">rank1</div>','<div class="image-label">rank2</div>'), $m);
    $sepn = 0;
    $m = preg_replace_callback('/<div class="data-separator"[^>]*>/', function($x) use (&$sepn){
        $sepn++;
        if ($sepn===1) return '<div class="data-separator" style="background-color:transparent">'; // data|rank1: no line (gap only)
        if ($sepn===2) return '<div class="data-separator">';                                       // rank2|mal1: single gray line
        return '<div class="data-separator" style="display:none">';                                 // #3,#4 around hidden mDOM/hDOM
    }, $m);
    // ---- paper layout: data label = length (kb), WHITE + right-aligned same as rank; straightened keeps raw+"(lin.)"; original hides Linearized ----
    // generic: applies to ANY folder that ships a _lengths.json (was hardcoded to basename=='Figure3')
    if (is_file("$base_path/_lengths.json")) {
        $isstr = (strpos($m, '>data(raw)</div>') !== false);  // straightened or not is read from the view itself (does a data(raw) panel exist), not from the
        // file name: the distribution strips the _straightened suffix, so the name no longer carries it
        if (!$isstr) {  // original: hide the Linearized panel (anchor on green label, BEFORE relabel strips color)
            $m = preg_replace('/<div class="comparison-item">(<div class="image-label" style="color:#27ae60")/', '<div class="comparison-item" style="display:none">$1', $m);
        }
        $LENJ = is_file("$base_path/_lengths.json") ? (json_decode(file_get_contents("$base_path/_lengths.json"), true) ?: array()) : array();
        if (isset($LENJ[$f]) && is_numeric($LENJ[$f])) {
            $lk = $LENJ[$f].' kb';
            // relabel to length; white (no inline color); SAME font as rank. flex justify-end pins the right
            // edge to box_right-15px (aligned with rank); a too-long label overflows LEFT into the black area (not clipped by the image on the right).
            $ds = 'style="display:flex;justify-content:flex-end;align-items:center"';
            $m = preg_replace('/<div class="image-label" style="color:#e67e22">data\(raw\)<\/div>/', '<div class="image-label" '.$ds.'><span style="white-space:nowrap">'.$lk.'</span></div>', $m);
            // the lin. label sits at the DNA contour height (~66%), not the image centre (50%): the
            // label-dots marker band pushes the molecule down
            $m = preg_replace('/<div class="image-label" style="color:#27ae60">Linearized<\/div>/', '<div class="image-label" style="position:relative;align-self:stretch"><span style="position:absolute;top:66%;right:15px;transform:translateY(-50%);white-space:nowrap">L '.$lk.'</span></div>', $m);
            $m = str_replace('<div class="image-label">data</div>', '<div class="image-label">'.$lk.'</div>', $m);
        }
    }
    // nav order from _summary.csv
    $bases = array(); $sum = "$base_path/_summary.csv";
    if (is_file($sum) && ($sh=fopen($sum,'r'))) { $hd=fgetcsv($sh); $mi=$hd?array_search('molecule',$hd):false;
        while(($rr=fgetcsv($sh))!==false){ if($mi!==false && !empty($rr[$mi]) &&
            (is_file("$base_path/{$rr[$mi]}_view_ytb.html") || is_file("$base_path/{$rr[$mi]}_view.html"))) $bases[]=$rr[$mi]; } fclose($sh); }
    if (!$bases){ foreach(array_merge(glob("$base_path/*_view_ytb.html"), glob("$base_path/*_view.html")) as $p)
                      $bases[]=preg_replace('/_view(_ytb)?\.html$/','',basename($p));
                  $bases=array_values(array_unique($bases)); sort($bases); }
    $i = array_search($f,$bases,true);
    $prev = ($i!==false && $i>0)?$bases[$i-1]:null; $next=($i!==false && $i!==count($bases)-1)?$bases[$i+1]:null;
    $pos = ($i===false)?1:$i+1; $tot=count($bases);
    $m = preg_replace('/(<span class="jump-num"[^>]*>)\d+(<\/span>\s*\/\s*)\d+/','${1}'.$pos.'${2}'.$tot,$m);
    $pct = $tot>0?round($pos*100.0/$tot,1):0;
    $m = preg_replace('/(progress-bar-fill"[^>]*width:)[\d.]+%/','${1}'.$pct.'%',$m);
    $m = preg_replace('/([A-Za-z0-9_.\-]+)_view_ytb\.html/','index.php?file=$1',$m);
    $m = str_replace('href="index.html"','href="index.php"',$m);
    $csvp = "$base_path/{$f}_matching_results.csv";
    if (is_file($csvp)) { $nt=hdom_table_from_csv($csvp);
        if ($nt) $m = preg_replace_callback('#<table[^>]*class="[^"]*matching-table[^"]*">.*?</table>#s', function($x)use($nt){return $nt;}, $m, 1); }
    $m = str_replace('<th>ref_name</th>','<th>chromosome</th>',$m);
    $m = str_replace('<th>log(P)</th>','<th>h_score</th>',$m);
    $mk = function($b){ return $b===null?'#':'index.php?file='.rawurlencode($b); };
    $m = preg_replace_callback('/(id="nav-prev"[^>]*href=")[^"]*(")/', function($x)use($prev,$mk){return $x[1].$mk($prev).$x[2];}, $m);
    $m = preg_replace_callback('/(id="nav-next"[^>]*href=")[^"]*(")/', function($x)use($next,$mk){return $x[1].$mk($next).$x[2];}, $m);
    // ref_map per-type show/hide toggle (local / mDOM / mal / hDOM) — persisted in localStorage
    $btn = '';
    foreach (array('local','mal') as $t) {
        $btn .= '<button type="button" data-rt="'.$t.'" data-off="0" onclick="tgR(this,\''.$t.'\')" '
             . 'style="margin:0 4px;padding:4px 12px;border:1px solid #888;border-radius:5px;'
             . 'cursor:pointer;background:#fff;font-weight:bold;">'.$t.'</button>';
    }
    $toggle = '<div style="text-align:center;padding:7px;background:#e9ecef;">'
        . '<span style="font-weight:bold;color:#555;margin-right:8px;">ref_map:</span>'.$btn.'</div>'
        . '<script>'
        . 'function _rmoff(){try{return JSON.parse(localStorage.getItem("refmap_off")||"{}")}catch(e){return{}}}'
        . 'function _rmapply(t,off){document.querySelectorAll(\'.comparison-item[data-rtype="\'+t+\'"]\').forEach(function(e){e.style.display=off?"none":"";});'
        . 'document.querySelectorAll(\'.matching-table tr[data-rtype="\'+t+\'"]\').forEach(function(e){e.style.display=off?"none":"";});}'
        . 'function tgR(b,t){var off=b.getAttribute("data-off")!=="1";_rmapply(t,off);'
        . 'b.setAttribute("data-off",off?"1":"0");b.style.opacity=off?"0.4":"1";'
        . 'var s=_rmoff();s[t]=off;localStorage.setItem("refmap_off",JSON.stringify(s));}'
        . 'document.addEventListener("DOMContentLoaded",function(){var s=_rmoff();'
        . 'document.querySelectorAll("button[data-rt]").forEach(function(b){var t=b.getAttribute("data-rt");'
        . 'if(s[t]){b.setAttribute("data-off","1");b.style.opacity="0.4";_rmapply(t,true);}});});'
        . '</script>';
    $m = preg_match('/<body[^>]*>/i', $m) ? preg_replace('/(<body[^>]*>)/i', '$1'.$toggle, $m, 1) : $toggle.$m;
    // navorder: if the listing was sorted, prev/next follow that order
    $navjs = '<script>(function(){try{var KEY="navorder:"+location.pathname.replace(/index\.php$/,"").replace(/\/$/,"");'
        . 'var order=JSON.parse(localStorage.getItem(KEY)||"null");if(!order||!order.length)return;'
        . 'var cur=new URLSearchParams(location.search).get("file");if(!cur)return;var i=order.indexOf(cur);if(i<0)return;'
        . 'function setn(id,b){var el=document.getElementById(id);if(!el)return;if(b==null){el.classList.add("disabled");el.removeAttribute("href");}else{el.classList.remove("disabled");el.setAttribute("href","index.php?file="+encodeURIComponent(b));}}'
        . 'function go(){setn("nav-prev",i>0?order[i-1]:null);setn("nav-next",i<order.length-1?order[i+1]:null);var jn=document.querySelector(".jump-num");if(jn)jn.textContent=(i+1);}'
        . 'if(document.getElementById("nav-next"))go();else document.addEventListener("DOMContentLoaded",go);}catch(e){}})();</script>';
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$navjs.'</body>',$m) : $m.$navjs;
    // ---- Cplus 'c'-key marking (C group only); state read server-side so it survives navigation ----
    if ($IS_C) {
        $marked = !empty(cplus_load($CPLUS_JSON)[$f]);
        $cjs = '<div id="cplus-banner" style="display:'.($marked?'block':'none').';position:fixed;top:0;left:0;right:0;'
            . 'z-index:99999;background:#b8860b;color:#fff;text-align:center;padding:6px;font-weight:bold;'
            . 'font-family:sans-serif;letter-spacing:.3px;">&#9873; marked &rarr; Cplus &nbsp;(press <b>c</b> to unmark)</div>'
            . '<script>(function(){var BASE='.json_encode($f).';var marked='.($marked?'true':'false').';'
            . 'var ban=document.getElementById("cplus-banner");function render(){if(ban)ban.style.display=marked?"block":"none";}'
            . 'document.addEventListener("keydown",function(e){if(e.key==="c"||e.key==="C"){'
            . 'var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():"";if(x==="input"||x==="textarea")return;'
            . 'e.preventDefault();var fd=new FormData();fd.append("action","toggle_cplus");fd.append("name",BASE);'
            . 'fetch("index.php",{method:"POST",body:fd}).then(function(r){return r.json();}).then(function(d){'
            . 'if(d&&d.ok){marked=d.marked;render();}}).catch(function(e){});}},false);render();})();</script>';
        $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$cjs.'</body>',$m) : $m.$cjs;
    }
    // ---- 't' trash marking (all groups); server-side (trash_marks.json), survives navigation ----
    $trashed = !empty(cplus_load($TRASH_JSON)[$f]);
    $tjs = '<div id="trash-banner" style="display:'.($trashed?'block':'none').';position:fixed;top:0;left:0;right:0;'
        . 'z-index:99998;background:#c0392b;color:#fff;text-align:center;padding:6px;font-weight:bold;'
        . 'font-family:sans-serif;letter-spacing:.3px;">&#128465; TRASH &mdash; marked for discard &nbsp;(press <b>t</b> to undo)</div>'
        . '<script>(function(){var BASE='.json_encode($f).';var tr='.($trashed?'true':'false').';'
        . 'var ban=document.getElementById("trash-banner");'
        . 'function render(){if(ban)ban.style.display=tr?"block":"none";'
        . 'document.documentElement.style.outline=tr?"5px solid #c0392b":"";document.documentElement.style.outlineOffset=tr?"-5px":"";}'
        . 'document.addEventListener("keydown",function(e){if(e.key==="t"||e.key==="T"){'
        . 'var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():"";if(x==="input"||x==="textarea")return;'
        . 'e.preventDefault();var fd=new FormData();fd.append("action","toggle_trash");fd.append("name",BASE);'
        . 'fetch("index.php",{method:"POST",body:fd}).then(function(r){return r.json();}).then(function(d){'
        . 'if(d&&d.ok){tr=d.trashed;render();}}).catch(function(e){});}},false);render();})();</script>';
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$tjs.'</body>',$m) : $m.$tjs;
    // ---- 'm' mutation/suspicious marking (all groups); server-side, survives navigation ----
    $mmarked = !empty(cplus_load($M_JSON)[$f]);
    $mjs = '<div id="m-banner" style="display:'.($mmarked?'block':'none').';position:fixed;top:0;left:0;right:0;'
        . 'z-index:99997;background:#6c3483;color:#fff;text-align:center;padding:6px;font-weight:bold;'
        . 'font-family:sans-serif;letter-spacing:.3px;">&#9888; M &mdash; mutation / suspicious (label diff) &nbsp;(press <b>m</b> to undo)</div>'
        . '<script>(function(){var BASE='.json_encode($f).';var mk='.($mmarked?'true':'false').';'
        . 'var ban=document.getElementById("m-banner");function render(){if(ban)ban.style.display=mk?"block":"none";}'
        . 'document.addEventListener("keydown",function(e){if(e.key==="m"||e.key==="M"){'
        . 'var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():"";if(x==="input"||x==="textarea")return;'
        . 'e.preventDefault();var fd=new FormData();fd.append("action","toggle_m");fd.append("name",BASE);'
        . 'fetch("index.php",{method:"POST",body:fd}).then(function(r){return r.json();}).then(function(d){'
        . 'if(d&&d.ok){mk=d.marked;render();}}).catch(function(e){});}},false);render();})();</script>';
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$mjs.'</body>',$m) : $m.$mjs;
    // ---- 'a' → A marking (all groups); server-side, survives navigation ----
    $amarked = !empty(cplus_load($A_JSON)[$f]);
    $ajs = '<div id="a-banner" style="display:'.($amarked?'block':'none').';position:fixed;top:0;left:0;right:0;'
        . 'z-index:99996;background:#1e8449;color:#fff;text-align:center;padding:6px;font-weight:bold;'
        . 'font-family:sans-serif;letter-spacing:.3px;">&#10004; marked &rarr; A group &nbsp;(press <b>a</b> to undo)</div>'
        . '<script>(function(){var BASE='.json_encode($f).';var mk='.($amarked?'true':'false').';'
        . 'var ban=document.getElementById("a-banner");function render(){if(ban)ban.style.display=mk?"block":"none";}'
        . 'document.addEventListener("keydown",function(e){if(e.key==="a"||e.key==="A"){'
        . 'var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():"";if(x==="input"||x==="textarea")return;'
        . 'e.preventDefault();var fd=new FormData();fd.append("action","toggle_a");fd.append("name",BASE);'
        . 'fetch("index.php",{method:"POST",body:fd}).then(function(r){return r.json();}).then(function(d){'
        . 'if(d&&d.ok){mk=d.marked;render();}}).catch(function(e){});}},false);render();})();</script>';
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$ajs.'</body>',$m) : $m.$ajs;
    // ---- 'k' → kink marking (all groups); server-side, survives navigation ----
    $kmarked = !empty(cplus_load($K_JSON)[$f]);
    $kjs = '<div id="k-banner" style="display:'.($kmarked?'block':'none').';position:fixed;top:0;left:0;right:0;'
        . 'z-index:99995;background:#8e44ad;color:#fff;text-align:center;padding:6px;font-weight:bold;'
        . 'font-family:sans-serif;letter-spacing:.3px;">&#12316; marked &rarr; kink (curvy) &nbsp;(press <b>k</b> to undo)</div>'
        . '<script>(function(){var BASE='.json_encode($f).';var mk='.($kmarked?'true':'false').';'
        . 'var ban=document.getElementById("k-banner");function render(){if(ban)ban.style.display=mk?"block":"none";}'
        . 'document.addEventListener("keydown",function(e){if(e.key==="k"||e.key==="K"){'
        . 'var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():"";if(x==="input"||x==="textarea")return;'
        . 'e.preventDefault();var fd=new FormData();fd.append("action","toggle_k");fd.append("name",BASE);'
        . 'fetch("index.php",{method:"POST",body:fd}).then(function(r){return r.json();}).then(function(d){'
        . 'if(d&&d.ok){mk=d.marked;render();}}).catch(function(e){});}},false);render();})();</script>';
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$kjs.'</body>',$m) : $m.$kjs;
    // ---- 'r' → AR marking (all groups); server-side (r_marks.json), survives navigation ----
    $rmarked = !empty(cplus_load($R_JSON)[$f]);
    $rjs = '<div id="r-banner" style="display:'.($rmarked?'block':'none').';position:fixed;top:0;left:0;right:0;'
        . 'z-index:99994;background:#d35400;color:#fff;text-align:center;padding:6px;font-weight:bold;'
        . 'font-family:sans-serif;letter-spacing:.3px;">&#9737; marked &rarr; AR (acrocentric p-arm / rDNA) &nbsp;(press <b>r</b> to undo)</div>'
        . '<script>(function(){var BASE='.json_encode($f).';var mk='.($rmarked?'true':'false').';'
        . 'var ban=document.getElementById("r-banner");function render(){if(ban)ban.style.display=mk?"block":"none";}'
        . 'document.addEventListener("keydown",function(e){if(e.key==="r"||e.key==="R"){'
        . 'var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():"";if(x==="input"||x==="textarea")return;'
        . 'e.preventDefault();var fd=new FormData();fd.append("action","toggle_r");fd.append("name",BASE);'
        . 'fetch("index.php",{method:"POST",body:fd}).then(function(r){return r.json();}).then(function(d){'
        . 'if(d&&d.ok){mk=d.marked;render();}}).catch(function(e){});}},false);render();})();</script>';

    // ---- image-quality banner, rendered from {base}_qc.json (written by hdom/qc.py) ----
    // Advisory only: the molecule is still aligned and still shown. The checks look at the IMAGE,
    // which is a different question from how the alignment scored, and one the score cannot answer.
    // Skipped when hDOM_view.py already baked a banner into the HTML, so the two never double up.
    $qcp = "$base_path/{$f}_qc.json";
    if (is_file($qcp) && strpos($m, 'data-qc-banner') === false) {
        $qc = json_decode((string)file_get_contents($qcp), true);
        if (is_array($qc) && empty($qc['pass'])) {
            $snr = isset($qc['green_snr']) ? number_format((float)$qc['green_snr'], 2) : '?';
            $med = (isset($qc['red_medroll']) && $qc['red_medroll'] !== null)
                 ? number_format((float)$qc['red_medroll'], 3) : 'untraceable';
            $gmin = $qc['green_snr_min'] ?? 2.5;
            $rmax = $qc['red_medroll_max'] ?? 0.33;
            $li = '';
            foreach (($qc['reasons'] ?? array()) as $r) $li .= '<li>'.htmlspecialchars($r).'</li>';
            $chip = function ($lab, $val, $bad) {
                return '<span style="display:inline-block;margin-right:14px;padding:2px 9px;border-radius:4px;'
                     . 'background:'.($bad ? '#c0392b' : '#7f8c8d').';color:#fff;font-weight:bold">'
                     . $lab.' '.$val.'</span>'; };
            $qcb = '<div data-qc-banner style="margin:10px 0;padding:12px 16px;border:1px solid #e0b400;'
                 . 'border-left:6px solid #e0b400;border-radius:5px;background:#fff8e1;color:#5a4500;'
                 . 'font:14px \'Segoe UI\',Tahoma,sans-serif;line-height:1.5">'
                 . '<div style="font-weight:bold;font-size:15px;margin-bottom:6px">'
                 . 'Image quality check failed: read this alignment with caution</div>'
                 . '<div style="margin-bottom:7px">'
                 . $chip('green SNR', $snr, !empty($qc['fail_green']))
                 . $chip('red medRoll', $med, !empty($qc['fail_red']))
                 . '<span style="color:#7a6a30;font-size:12.5px">thresholds: green SNR &ge; '.$gmin
                 . ', red medRoll &le; '.$rmax.'</span></div>'
                 . '<ul style="margin:0 0 7px 20px;padding:0">'.$li.'</ul>'
                 . '<div style="font-size:12.5px;color:#7a6a30">These checks look at the image only, not at '
                 . 'the alignment. The result below was still computed and is shown unchanged; the alignment '
                 . 'score cannot tell you that the labels it ranked on were mostly background.</div></div>';
            if (preg_match('/<div class="eval-toolbar"/', $m)) {
                $m = preg_replace('/(<div class="eval-toolbar")/', $qcb.'$1', $m, 1);
            } elseif (preg_match('/<body[^>]*>/i', $m)) {
                $m = preg_replace('/(<body[^>]*>)/i', '$1'.$qcb, $m, 1);
            } else { $m = $qcb.$m; }
        }
    }
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$rjs.'</body>',$m) : $m.$rjs;
    echo $m; exit;
}

// ---- listing (A_group style: sortable + navorder persistence) ----
$rows = array();
if (is_file("$base_path/_summary.csv") && ($fp=fopen("$base_path/_summary.csv","r"))) {
    $hd=fgetcsv($fp); while(($r=fgetcsv($fp))!==false) $rows[]=array_combine($hd,array_pad($r,count($hd),'')); fclose($fp); }
$g = basename($base_path);
$LEN = is_file("$base_path/_lengths.json") ? (json_decode(file_get_contents("$base_path/_lengths.json"),true) ?: array()) : array();
echo "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1.0'>"
   . "<title>$g — local-align comparison</title><style>"
   . "body{font-family:'Segoe UI',Tahoma,sans-serif;margin:20px;background:#f8f9fa;color:#333}"
   . "h1{color:#0056b3;border-bottom:2px solid #0056b3;padding-bottom:10px}"
   . "table{width:100%;border-collapse:collapse;background:#fff;box-shadow:0 1px 3px rgba(0,0,0,.2);margin-top:16px}"
   . "th,td{padding:8px 12px;text-align:left;border-bottom:1px solid #ddd;font-size:.92em}"
   . "th{background:#0056b3;color:#fff;cursor:pointer;white-space:nowrap}th:hover{background:#004494}"
   . "th::after{content:' \\2195';font-size:.8em;color:#ccc}tr:hover{background:#f1f1f1}"
   . "a{color:#0056b3;text-decoration:none;font-weight:bold}a:hover{text-decoration:underline}"
   . ".small{font-weight:normal;color:#888}.toolbar a{display:inline-block;padding:4px 12px;background:#305496;color:#fff;border-radius:4px;font-size:.9em;margin-right:6px}</style></head><body>";
echo "<h1>$g <span class='small'>(".count($rows)." molecules) &middot; local-align 4-method comparison</span></h1>";
$cmarks = $IS_C ? cplus_load($CPLUS_JSON) : array();
echo "<div class='toolbar'><a href='../'>&larr; DOM_local</a><a href='_summary.csv'>_summary.csv</a>";
if ($IS_C) echo "<span style='display:inline-block;padding:4px 12px;border-radius:4px;background:#b8860b;color:#fff;"
    . "font-weight:bold;font-size:.9em'>&#9873; <span id='mvn'>".count($cmarks)."</span> marked &rarr; Cplus</span>"
    . "<span class='small' style='margin-left:8px'>(open a molecule, press <b>c</b> to mark; run move_cplus.py to move)</span>";
echo "</div>";
echo "<table id='t'><thead><tr>";
$scoremethod = 'local';  // ML_label tree: _summary score = local rank1 logP (local_cc_r x F1)  // local tree -> local logP
$othermethod = ($scoremethod==='local') ? 'mDOM' : 'local';  // other method's classification column
$hasmc = false;  // mDOM class column removed (per request)
$cols  = array('#','Molecule','len(kb)','H_score ('.$scoremethod.' logP)','local rank1','mal rank1');
$types = array('s','s','n','n','s','s','s','s','s');
if ($hasmc){ $cols[]=$othermethod.' class'; $types[]='s'; }
$hashom = !empty($rows) && array_key_exists('hom12',$rows[0]);
if ($hashom){ $cols[]='homolog 1,2'; $types[]='s'; }
$haschg = !empty($rows) && array_key_exists('change',$rows[0]);
if ($haschg){ $cols[]='change'; $types[]='s'; }
$hasog = !empty($rows) && array_key_exists('localgrp',$rows[0]);
if ($hasog){ $cols[]='origin (local)'; $types[]='s'; }
if ($IS_C){ $cols[]='Cplus'; $types[]='s'; }
$cols[]='view'; $types[]='';
foreach($cols as $k=>$c){ echo ($c==='view') ? "<th>view</th>" : "<th onclick=\"sortTable($k,'{$types[$k]}')\">".htmlspecialchars($c)."</th>"; }
echo "</tr></thead><tbody>";
$i=0; foreach($rows as $r){ $i++; $mol=$r['molecule']; $v='index.php?file='.rawurlencode($mol);
 echo "<tr data-mol='".htmlspecialchars($mol)."'><td>$i</td>"
    ."<td><a href='$v'>".htmlspecialchars($mol)."</a></td>"
    ."<td style='text-align:right'>".(isset($LEN[$mol])?$LEN[$mol]:'')."</td>"
    ."<td>".htmlspecialchars($r['score']??'')."</td>"
    ."<td>".htmlspecialchars($r['local1']??'')."</td>"
    ."<td>".htmlspecialchars($r['mal1']??'')."</td>"
    
    .($hasmc ? "<td>".htmlspecialchars($r['mdom_class']??'')."</td>" : "")
    .($hashom ? "<td style='text-align:center;font-weight:bold;color:".((($r['hom12']??'')==='Y')?'#27ae60':'#bbb')."'>".((($r['hom12']??'')==='Y')?'&#10004;':'&ndash;')."</td>" : "")
    .($haschg ? "<td style='text-align:center;font-weight:bold;color:".(($r['change']??'')==='GAIN'?'#27ae60':(($r['change']??'')==='BREAK'?'#e74c3c':'#bbb'))."'>".htmlspecialchars($r['change']??'')."</td>" : "")
    .($hasog ? "<td style='text-align:center;font-weight:bold'>".htmlspecialchars($r['localgrp']??'')."</td>" : "")
    .($IS_C ? "<td class='cpm' data-mol='".htmlspecialchars($mol)."' style='text-align:center;color:#b8860b;font-weight:bold'>".(!empty($cmarks[$mol])?'&#9873;':'')."</td>" : "")
    ."<td><a href='$v'>view</a></td></tr>";
}
echo "</tbody></table><script>"
  . "function sortTable(n,type){var t=document.getElementById('t'),tb=t.tBodies[0];var rows=Array.from(tb.rows);"
  . "var th=t.tHead.rows[0].cells[n];var dir=th.getAttribute('data-dir')==='asc'?'desc':'asc';"
  . "Array.from(t.tHead.rows[0].cells).forEach(function(c){if(c!==th)c.removeAttribute('data-dir');});th.setAttribute('data-dir',dir);"
  . "rows.sort(function(a,b){var x=a.cells[n].textContent.trim(),y=b.cells[n].textContent.trim();"
  . "if(type==='n'){x=parseFloat(x);y=parseFloat(y);if(isNaN(x))x=Infinity;if(isNaN(y))y=Infinity;}"
  . "if(x<y)return dir==='asc'?-1:1;if(x>y)return dir==='asc'?1:-1;return 0;});"
  . "rows.forEach(function(r){tb.appendChild(r);});saveNavOrder();}"
  . "function navKey(){return 'navorder:'+location.pathname.replace(/index\\.php$/,'').replace(/\\/$/,'');}"
  . "function saveNavOrder(){try{var tb=document.getElementById('t').tBodies[0];"
  . "var order=Array.from(tb.rows).map(function(r){return r.getAttribute('data-mol');}).filter(Boolean);"
  . "localStorage.setItem(navKey(),JSON.stringify(order));}catch(e){}}saveNavOrder();"
  . "</script></body></html>";
