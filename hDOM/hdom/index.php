<?php
register_shutdown_function(function(){ if(($_SERVER['REQUEST_METHOD']??'')==='POST' || ($_GET['file']??'')==='') return; $h=$GLOBALS['HDOM_UP_HREF'] ?? 'index.php'; echo '<a href="'.$h.'" style="position:fixed;top:6px;left:6px;z-index:2147483647;background:#34495e;color:#fff;padding:6px 14px;border-radius:5px;text-decoration:none;font:bold 13px sans-serif;box-shadow:0 1px 4px rgba(0,0,0,.4)">&#8593; Up</a>'; });
// 4-method local-align comparison index. Serves <base>_view_ytb.html, or <base>_view.html:
// the analysis server writes the first, run_pipeline.py writes the second, and one viewer has to
// read both or a released run lands on a 404.
// and re-renders the matching table (local / mDOM / mal / pcc + H) from the CSV.
$base_path = __DIR__;
// Helpers that point prev/next, Index and Up at the Total list when the page was entered from it (?from=total).
foreach (array(__DIR__.'/../_total_nav.php', __DIR__.'/../../_total_nav.php') as $__tn)
    if (is_file($__tn)) { require_once $__tn; break; }
foreach (array(__DIR__.'/../_summary_sync.php', __DIR__.'/../../_summary_sync.php') as $__ss)
    if (is_file($__ss)) { require_once $__ss; break; }
header('Cache-Control: no-store, no-cache, must-revalidate');

// ---- mark-root: all *_marks.json live at the DOM_local root (the ancestor dir named 'local'),
// NOT dirname($base_path). Nested folders (zave/fresh_reclass/A_to_C, zave/C_total/*, ...) are
// several levels below 'local'; without this their marks would try to write into their own parent
// (non-existent, apache-unwritable) json and silently vanish, and reclassify (which reads local/*)
// would never see them. Walk up to 'local' so every folder shares the same canonical mark files.
$__p = $base_path;
while ($__p && $__p !== '/' && $__p !== '.' && basename($__p) !== 'local' && basename($__p) !== 'ML_label_local_ccr' && basename($__p) !== 'hDOM_HG002_v1.2_1bp') { $__p = dirname($__p); }
$MARK_ROOT = (basename($__p) === 'local' || basename($__p) === 'ML_label_local_ccr' || basename($__p) === 'hDOM_HG002_v1.2_1bp') ? $__p : dirname($base_path);

// ---- Cplus marking (C group only): server-side flag, survives navigation ----
// Marks live at local/cplus_marks.json (root, survives reclassify wipe) and must be
// pre-created world-writable (666) by reclassify so apache can overwrite them.
// C-marking enabled for the C group, C_total subfolders, and ANY fresh_reclass transition bucket
// (buckets are per-transition like C_to_B / A_to_C, so match by parent dir rather than listing each).
$IS_C = in_array(basename($base_path), array('C', 'C_longmolecules', 'C_longer_than_200kb', 'candidate_recovered', 'linearization', 'others', 'Long_Molecules', 'C_freshscan_income'), true)
     || basename(dirname($base_path)) === 'fresh_reclass';
$CPLUS_JSON = $MARK_ROOT . '/cplus_marks.json';
// When the Total list was sorted, the neighbours follow that order. The Total page leaves
// navorder:total = {m:[molecule...], p:{molecule:"group/sub/"}} in localStorage; the paths are
// rebuilt against this page's relative depth (rel). A stored list whose length no longer matches
// the current total is stale and is ignored.
function hdom_total_navjs($T){
    return '<script>(function(){try{var N='.(int)$T['n'].';var REL='.json_encode($T['rel']).';'
      . 'var o=JSON.parse(localStorage.getItem("navorder:total")||"null");'
      . 'if(!o||!o.m||o.m.length!==N||!o.p)return;'
      . 'var cur=new URLSearchParams(location.search).get("file");if(!cur)return;'
      . 'var i=o.m.indexOf(cur);if(i<0)return;'
      . 'function href(b){return REL+o.p[b]+"index.php?file="+encodeURIComponent(b)+"&from=total";}'
      . 'function setn(id,b){var el=document.getElementById(id);if(!el)return;'
      . 'if(b==null){el.classList.add("disabled");el.removeAttribute("href");}'
      . 'else{el.classList.remove("disabled");el.setAttribute("href",href(b));}}'
      . 'function go(){setn("nav-prev",i>0?o.m[i-1]:null);setn("nav-next",i<N-1?o.m[i+1]:null);'
      . 'var jn=document.querySelector(".jump-num");if(jn)jn.textContent=(i+1);}'
      . 'if(document.getElementById("nav-next"))go();else document.addEventListener("DOMContentLoaded",go);'
      . '}catch(e){}})();<\/script>';
}
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
    // shift_gap is the mean disagreement in px between the best shift of neighbouring windows.
    // w_coh = 1 - shift_gap/(2*MAXS) is the factor the score is multiplied by, 1 being perfectly
    // smooth. The table carries the factor only.
    $scf = function($r){
        $w = $r['w_coh'] ?? ''; 
        return ($w === '' || !is_numeric($w)) ? '' : sprintf('%.2f', (float)$w); };

    // Row colours: hDOM = method-local (#e8f4f8), mal = method-mal (#fff4e6). The mal header is a quieter tan.
    $css = '<style>'
      .'.mm{font-size:1.10em}'                      // matches the size of the saved figure
      .'.mm .hdr-hdom td{background:#3498db;color:#fff;font-weight:bold}'
      .'.mm .hdr-mal td{background:#a67c52;color:#fff;font-weight:bold}'
      .'.mm .osep td{background:#34495e;color:#fff;font-weight:bold;text-align:left}'
      .'.mm .hdr-hdom,.mm .hdr-mal,.mm .osep{cursor:default}'
      .'.mm td:nth-child(n+6),.mm th:nth-child(n+6){width:7.30%;min-width:56px;text-align:center;box-sizing:border-box}'
      .'.mm td:nth-child(12),.mm th:nth-child(12){width:8.0%;min-width:64px}'  // scale..h_score and mal share one uniform width (7 columns x 9.5% = 66.5%)
      .'.mm thead th,.mm .hdr-hdom td,.mm .hdr-mal td{padding-top:4px !important;padding-bottom:4px !important;line-height:1.1 !important}'  // header row: same font size, thinner vertically
      // Data rows (hDOM/mal) only, at 0.9x height: padding 8.1+8.1 plus line 17.3 = 33.5px.
      // The header row (thead th), the group headers (hdr-*) and the Other ranks separator (osep)
      // are excluded so the rules above still hold for them.
      .'.mm tbody tr:not(.hdr-hdom):not(.hdr-mal):not(.osep) td{padding-top:7px !important;padding-bottom:7px !important;line-height:1.42 !important}'
      .'</style>';

    // hDOM has 12 columns (cc_rg and misc dropped). mal leaves the cc_r slot (7th) empty and lines matched..m_score up.
    $HDR_H = array('rank','chromosome','method','aligned_bp','dir','scale','cc_r','cc_g','cc_rg2','cc_rw','F1_g','coh_w','h_score');
    $HDR_M = array('rank','chromosome','method','aligned_bp','dir','scale','label','matched','FP','FN','F1','sizing','m_score');
    $hrow = function($cols,$cls){ $t='    <tr class="'.$cls.'">'; foreach($cols as $c) $t.='<td>'.htmlspecialchars($c).'</td>'; return $t."</tr>\n"; };
    $rLoc = function($pair) use($f2c,$fec,$bpf,$drf,$scf){
        list($i,$r)=$pair; $rk=(string)($r['rank']??'');
        $cl='method-local'.((is_numeric($rk)&&((int)$rk===1||(int)$rk===2))?' rank-top':'');
        $tds='<td>rank'.htmlspecialchars(is_numeric($rk)?(string)(int)$rk:$rk).'</td>'
            .'<td>'.htmlspecialchars((string)($r['ref_name']??'')).'</td><td>hDOM</td>'
            .'<td>'.$bpf($r).'</td><td>'.$drf($r).'</td><td>'.$f2c($r,'scale').'</td>'
            .'<td>'.$f2c($r,'cc_r').'</td><td>'.$f2c($r,'cc_g').'</td><td>'.$fec($r,'cc_rg2').'</td>'
            .'<td>'.$f2c($r,'local_cc_r').'</td><td>'.$f2c($r,'F1').'</td>'
            .'<td>'.$scf($r).'</td><td>'.$f2c($r,'log(P)').'</td>';
        return '    <tr class="'.$cl.'" data-rtype="local" data-sel="'.($i+1).'">'.$tds."</tr>\n";
    };
    $rMal = function($pair) use($f2c,$txf,$bpf,$drf){
        list($i,$r)=$pair; $rk=(string)($r['rank']??'');
        $cl='method-mal'.((is_numeric($rk)&&((int)$rk===1||(int)$rk===2))?' rank-top':'');
        // maligner's mal_matched counts matched chunks; K chunks join K-1 labels. matched, F1
        // and label are all put on a per-label basis here (FP = query_misses and FN = ref_misses
        // are label counts already). This is the TP definition Total/metrics.json uses.
        $mChunk=(int)($r['mal_matched']??0); $mM=max(0,$mChunk-1);
        $mFP=(int)($r['mal_FP']??0); $mFN=(int)($r['mal_FN']??0);
        $mden=2*$mM+$mFP+$mFN; $mf1=($mden>0)?number_format(2*$mM/$mden,2):'';   // F1 = 2M/(2M+FP+FN)
        $mlab=(($r['mal_matched']??'')!=='' && ($r['mal_FP']??'')!=='') ? (string)($mM+$mFP) : '';
        $tds='<td>m'.htmlspecialchars(is_numeric($rk)?(string)(int)$rk:$rk).'</td>'
            .'<td>'.htmlspecialchars((string)($r['ref_name']??'')).'</td><td>mal</td>'
            .'<td>'.$bpf($r).'</td><td>'.$drf($r).'</td><td>'.$f2c($r,'scale').'</td>'
            // label = the labels this molecule entered maligner with = matched + query_misses (FP)
            .'<td>'.$mlab.'</td>'
            .'<td>'.(($r['mal_matched']??'')!=='' ? (string)$mM : '').'</td>'
            .'<td>'.$txf($r,'mal_FP').'</td><td>'.$txf($r,'mal_FN').'</td>'
            .'<td>'.$mf1.'</td><td>'.$f2c($r,'mal_sizing').'</td><td>'.$f2c($r,'m_score').'</td>';
        return '    <tr class="'.$cl.'" data-rtype="mal" data-sel="'.($i+1).'">'.$tds."</tr>\n";
    };

    // One table, 12 columns laid out automatically. thead is the hDOM header.
    $t = '<table border="0" class="dataframe matching-table mm">'."\n  <thead>\n    <tr>\n";
    foreach($HDR_H as $c) $t.='      <th>'.htmlspecialchars($c)."</th>\n"; $t.="    </tr>\n  </thead>\n  <tbody>\n";
    foreach (array_slice($loc,0,3) as $p) $t.=$rLoc($p);          // top: hDOM top-3
    $t.=$hrow($HDR_M,'hdr-mal');
    foreach (array_slice($mal,0,3) as $p) $t.=$rMal($p);          //          mal top-3
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
    // Homolog panel labels: only the long ones that carry a chromosome name (H_local chr05_p and
    // the like) are shrunk. Taking a short label such as 'H_local' down to 9px only makes its font
    // differ from the neighbouring panels.
    $m = preg_replace_callback('/<div class="image-label">(H_[^<]*)<\/div>/', function($x){
        return (strpos($x[1], 'chr') !== false)
            ? '<div class="image-label" style="font-size:9px;padding:0 4px;letter-spacing:-.2px">'.$x[1].'</div>'
            : $x[0]; }, $m);
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
    // The rule that hides the mDOM/hDOM separators also catches the one before the homolog panel,
    // which leaves the H panel flush against the panel above. That single position keeps the same
    // grey line rank2|mal1 has.
    $m = preg_replace('/<div class="data-separator" style="display:none"><\/div>(<div class="comparison-item"[^>]*><div class="image-label"[^>]*>H_)/',
                      '<div class="data-separator"></div>$1', $m);
    // ---- paper layout: data label = length (kb), WHITE + right-aligned same as rank; straightened keeps raw+"(lin.)"; original hides Linearized ----
    // generic: applies to ANY folder that ships a _lengths.json (was hardcoded to basename=='Figure3')
    if (is_file("$base_path/_lengths.json")) {
        $isstr = (strpos($m, '>data(raw)</div>') !== false);  // straightened is read from the view's content (a data(raw) panel is present), not from the file name, which carries no _straightened suffix
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
            // the lin. label sits at the DNA contour height (~66%), not the image centre: the label-dots marker strip pushes the molecule down
            $m = preg_replace('/<div class="image-label" style="color:#27ae60">Linearized<\/div>/', '<div class="image-label" style="position:relative;align-self:stretch"><span style="position:absolute;top:66%;right:15px;transform:translateY(-50%);white-space:nowrap">L '.$lk.'</span></div>', $m);
            $m = str_replace('<div class="image-label">data</div>', '<div class="image-label">'.$lk.'</div>', $m);
        }
    }
    // Entered from the Total list: that list defines the neighbours, the counter and Up.
    $TCTX = function_exists('hdom_total_ctx') ? hdom_total_ctx($base_path,$f) : null;
    if ($TCTX && $TCTX['i']!==false) $GLOBALS['HDOM_UP_HREF'] = $TCTX['up'];
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
    if ($TCTX && $TCTX['i']!==false) { $pos=$TCTX['pos']; $tot=$TCTX['n']; }
    $m = preg_replace('/(<span class="jump-num"[^>]*>)\d+(<\/span>\s*\/\s*)\d+/','${1}'.$pos.'${2}'.$tot,$m);
    $pct = $tot>0?round($pos*100.0/$tot,1):0;
    $m = preg_replace('/(progress-bar-fill"[^>]*width:)[\d.]+%/','${1}'.$pct.'%',$m);
    $m = preg_replace('/([A-Za-z0-9_.\-]+)_view_ytb\.html/','index.php?file=$1',$m);
    $m = str_replace('href="index.html"','href="index.php"',$m);
    if ($TCTX) $m = str_replace('class="index-link" href="index.php"','class="index-link" href="'.htmlspecialchars($TCTX['up'],ENT_QUOTES,'UTF-8').'"',$m);
    $csvp = "$base_path/{$f}_matching_results.csv";
    if (is_file($csvp)) { $nt=hdom_table_from_csv($csvp);
        if ($nt) $m = preg_replace_callback('#<table[^>]*class="[^"]*matching-table[^"]*">.*?</table>#s', function($x)use($nt){return $nt;}, $m, 1); }
    $m = str_replace('<th>ref_name</th>','<th>chromosome</th>',$m);
    $m = str_replace('<th>log(P)</th>','<th>h_score</th>',$m);
    // --- Select rank: drop mDOM-era d* keys, add rank steppers ---
    $m = preg_replace('/<option value="d\d+">\s*/', '', $m);
    $m = preg_replace_callback('/(LABEL_MAP\s*=\s*\{)([^}]*)(\})/', function($x){
        return $x[1].rtrim(rtrim(preg_replace('/"d\d+"\s*:\s*\d+\s*,?\s*/', '', $x[2])), ',').$x[3]; }, $m, 1);
    $m = str_replace('placeholder="e.g. d1, m3, p2"', 'placeholder="e.g. l1, m3, p2"', $m);
    $m = str_replace('Enter d1.., m1.., p1.. to display.', 'Enter l1.., m1.., p1.. to display.', $m);
    $sbtn = 'line-height:1;padding:1px 7px;font-size:9px;border:1px solid #888;background:#f2f2f2;cursor:pointer';
    $sstep = '<span id="sel-step" style="display:inline-flex;flex-direction:column;margin-left:5px;vertical-align:middle">'
           . '<button type="button" onclick="selStep(-1)" title="prev rank (Up arrow)" style="'.$sbtn.';border-bottom:none;border-radius:4px 4px 0 0">&#9650;</button>'
           . '<button type="button" onclick="selStep(1)" title="next rank (Down arrow)" style="'.$sbtn.';border-radius:0 0 4px 4px">&#9660;</button></span>';
    $m = preg_replace('/(<input type="text" id="sel-input"[^>]*>)/', '$1'.$sstep, $m, 1);
    $sjs = '<script>function selStep(d){var i=document.getElementById("sel-input");if(!i)return;'
         . 'var L=Array.prototype.map.call(document.querySelectorAll("#sel-rank-list option"),function(o){return o.value;});'
         . 'if(!L.length)return;var v=(i.value||"").trim().toLowerCase();var t=v.match(/^([a-z]+)(\d+)$/);'
         . 'var p=t?t[1]:"l";var s=L.filter(function(x){var q=x.match(/^([a-z]+)(\d+)$/);return q&&q[1]===p;});'
         . 'if(!s.length)s=L;var j=s.indexOf(v);j=(j<0)?0:j+d;if(j<0)j=0;if(j>=s.length)j=s.length-1;'
         . 'i.value=s[j];i.dispatchEvent(new Event("input",{bubbles:true}));i.dispatchEvent(new Event("change",{bubbles:true}));}'
         . 'document.addEventListener("keydown",function(e){if(e.target&&e.target.id==="sel-input"){'
         . 'if(e.key==="ArrowUp"){e.preventDefault();selStep(-1);}else if(e.key==="ArrowDown"){e.preventDefault();selStep(1);}}},true);</script>';
    $m = (strpos($m,'</body>')!==false) ? str_replace('</body>',$sjs.'</body>',$m) : $m.$sjs;
    $mk = function($b){ return $b===null?'#':'index.php?file='.rawurlencode($b); };
    $purl = $TCTX ? ($TCTX['prev']===null?'#':$TCTX['prev']) : $mk($prev);
    $nurl = $TCTX ? ($TCTX['next']===null?'#':$TCTX['next']) : $mk($next);
    if ($TCTX) {
        // In Total order the disabled state and href baked in per folder can both be wrong: at a
        // group boundary a molecule is the folder's last but not Total's. The tags are rewritten
        // whole.
        foreach (array(array('nav-prev',$TCTX['prev'],'Previous molecule (Left arrow)','&#9664;'),
                       array('nav-next',$TCTX['next'],'Next molecule (Right arrow)','&#9654;')) as $nb) {
            list($id,$u,$ti,$gl) = $nb;
            $tag = '<a id="'.$id.'" class="nav-btn'.($u===null?' disabled':'').'"'
                 . ($u===null?'':' href="'.htmlspecialchars($u,ENT_QUOTES,'UTF-8').'"')
                 . ' title="'.$ti.'">'.$gl.'</a>';
            $m = preg_replace('#<a id="'.$id.'"[^>]*>.*?</a>#s', $tag, $m, 1);
        }
    } else {
        $m = preg_replace_callback('/(id="nav-prev"[^>]*href=")[^"]*(")/', function($x)use($purl){return $x[1].htmlspecialchars($purl,ENT_QUOTES,'UTF-8').$x[2];}, $m);
        $m = preg_replace_callback('/(id="nav-next"[^>]*href=")[^"]*(")/', function($x)use($nurl){return $x[1].htmlspecialchars($nurl,ENT_QUOTES,'UTF-8').$x[2];}, $m);
    }
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
    $navjs = $TCTX ? hdom_total_navjs($TCTX) : '<script>(function(){try{var TOTN='.count($bases).';var KEY="navorder:"+location.pathname.replace(/index\.php$/,"").replace(/\/$/,"");'
        . 'var order=JSON.parse(localStorage.getItem(KEY)||"null");if(!order||!order.length)return;'
        // A stored order from a listing with a different molecule count is stale (files added or
        // removed since); drop it and keep the server-rendered prev/next rather than navigating
        // by a list that no longer describes this folder.
        . 'if(order.length!==TOTN){localStorage.removeItem(KEY);return;}'
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

    // ---- recovery banner, rendered from {base}_recovery.json ----
    // Same slot and same shape as the image-quality banner below, in blue rather than amber: this one
    // is not a warning about the molecule, it states that the alignment shown was produced under a
    // configuration other than the deployed default, and which one. Molecules without the sidecar are
    // ordinary deployed results and get no banner.
    $rcp = "$base_path/{$f}_recovery.json";
    if (is_file($rcp) && strpos($m, 'data-recovery-banner') === false) {
        $rc = json_decode((string)file_get_contents($rcp), true);
        if (is_array($rc) && !empty($rc['title'])) {
            $rchip = function ($lab, $val) {
                return '<span style="display:inline-block;margin-right:14px;padding:2px 9px;border-radius:4px;'
                     . 'background:#3d6fa5;color:#fff;font-weight:bold">'
                     . htmlspecialchars($lab).' '.htmlspecialchars($val).'</span>'; };
            $ch = '';
            foreach (($rc['chips'] ?? array()) as $cpair) {
                if (is_array($cpair) && count($cpair) >= 2) $ch .= $rchip($cpair[0], $cpair[1]); }
            $na = (int)($rc['n_add'] ?? 0); $nr = (int)($rc['n_rem'] ?? 0);
            if ($na || $nr) $ch .= $rchip('label edits', sprintf('+%d / -%d', $na, $nr));
            $og = $rc['origin'] ?? '';
            $rcb = '<div data-recovery-banner style="margin:10px 0;padding:12px 16px;border:1px solid #3d6fa5;'
                 . 'border-left:6px solid #3d6fa5;border-radius:5px;background:#eef4fb;color:#22405f;'
                 . 'font:14px \'Segoe UI\',Tahoma,sans-serif;line-height:1.5">'
                 . '<div style="font-weight:bold;font-size:15px;margin-bottom:6px">'
                 . htmlspecialchars($rc['title']).' &mdash; not the deployed configuration</div>'
                 . ($ch !== '' ? '<div style="margin-bottom:7px">'.$ch.'</div>' : '')
                 . '<div style="margin-bottom:7px">'.htmlspecialchars($rc['body'] ?? '').'</div>'
                 . '<div style="font-size:12.5px;color:#4a6a8c">The deployed run leaves this molecule without a '
                 . 'homolog pair'.($og !== '' ? ', in '.htmlspecialchars($og) : '').'. h_score comes from the same '
                 . 'K=5000 null as every other molecule, but the candidate list it was ranked against is not the '
                 . 'default one, so read it against this molecule&rsquo;s own alternatives rather than across molecules.'
                 . '</div></div>';
            if (preg_match('/<div class="eval-toolbar"/', $m)) {
                $m = preg_replace('/(<div class="eval-toolbar")/', $rcb.'$1', $m, 1);
            } elseif (preg_match('/<body[^>]*>/i', $m)) {
                $m = preg_replace('/(<body[^>]*>)/i', '$1'.$rcb, $m, 1);
            } else { $m = $rcb.$m; }
        }
    }

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
    // B2 header layout: the ref_map bar, the nav arrows and the Class toolbar are rearranged and
    // the duplicated subtitle is dropped. The profile and cc graphs and the h1 analysis table are
    // not taken.
    $__bl = __DIR__;
    for ($__i = 0; $__i < 4 && !is_file($__bl.'/_b2_layout.php'); $__i++) $__bl = dirname($__bl);
    if (is_file($__bl.'/_b2_layout.php')) {
        require_once $__bl.'/_b2_layout.php';
        $__lay = hdom_b2_layout();
        $m = (strpos($m,'</body>') !== false) ? str_replace('</body>', $__lay.'</body>', $m) : $m.$__lay;
    }
    echo $m; exit;
}

// ---- listing (hDOM Database card style: sortable + navorder persistence) ----
// When a molecule's CSV is newer than the summary, after a re-rank or an h_score recomputation,
// score and rank1 are rebuilt on the spot. _summary.csv is updated too where the process can
// write; otherwise the fresh values are used for the page only.
$rows = array();
if (function_exists('hdom_summary_rows')) $rows = hdom_summary_rows($base_path) ?: array();
if (!$rows && is_file("$base_path/_summary.csv") && ($fp=fopen("$base_path/_summary.csv","r"))) {
    $hd=fgetcsv($fp); while(($r=fgetcsv($fp))!==false) $rows[]=array_combine($hd,array_pad($r,count($hd),'')); fclose($fp); }
$g = basename($base_path);
// Title, description and colour come from the landing page's $GROUPS block, so a folder name never
// surfaces as a title. A sub-set (C_total/Information_rich ...) is described by its parent
// container's $SUBS entries and takes the parent group's colour.
$GMETA = array(); $SUBMETA = array();
for ($up=$base_path, $i=0; $i<4; $i++) {
    $cand = dirname($up).'/index.php';
    if (is_file($cand)) {
        $src = file_get_contents($cand);
        if (!$SUBMETA && preg_match_all('/array\(\'([A-Za-z0-9_]+)\',\'([^\']*)\',\s*\'([^\']*)\'\)/', $src, $sm, PREG_SET_ORDER))
            foreach ($sm as $s) $SUBMETA[$s[1]] = array($s[2], $s[3]);
        if (preg_match('/\$GROUPS *= *json_decode\(<<<\'JSON\'(.*?)JSON/s', $src, $mm)) {
            $j = json_decode(trim($mm[1]), true);
            if (is_array($j)) { $GMETA = $j; break; }
        }
    }
    $up = dirname($up);
}
$pg = basename(dirname($base_path));
$gm = (isset($GMETA[$g]) && is_array($GMETA[$g])) ? $GMETA[$g] : null;
$pm = (isset($GMETA[$pg]) && is_array($GMETA[$pg])) ? $GMETA[$pg] : null;
$disp  = $gm ? $gm[1] : (isset($SUBMETA[$g]) ? $SUBMETA[$g][0] : str_replace('_', ' ', $g));
$gdesc = $gm ? (string)($gm[2] ?? '') : (isset($SUBMETA[$g]) ? $SUBMETA[$g][1] : '');
$gcol  = $gm ? (string)($gm[4] ?? '') : ($pm ? (string)($pm[4] ?? '') : '');
$back  = $pm ? $pm[1] : 'hDOM Database';
$LEN = is_file("$base_path/_lengths.json") ? (json_decode(file_get_contents("$base_path/_lengths.json"),true) ?: array()) : array();
// Label-only folders (SGM and its Figure 2 subset) have no hDOM result: local rank1 is empty on every
// row and the _summary score column is the maligner m_score.
$MALONLY = !empty($rows);
foreach ($rows as $r) if (($r['local1'] ?? '') !== '') { $MALONLY = false; break; }
$esc = function($s){ return htmlspecialchars((string)$s, ENT_QUOTES, 'UTF-8'); };
$cols = array(array('#','n','idx'), array('Molecule','s','mol'), array('len (kb)','n','num'),
              $MALONLY ? array('m_score','n','num') : array('h_score (logP)','n','num'));
if (!$MALONLY) $cols[] = array('hDOM rank1','s','');
$cols[] = array('mal rank1','s','');
$hashom = !empty($rows) && array_key_exists('hom12',$rows[0]);    if ($hashom) $cols[] = array('homolog 1,2','s','ctr');
$haschg = !empty($rows) && array_key_exists('change',$rows[0]);   if ($haschg) $cols[] = array('change','s','ctr');
$hasog  = !empty($rows) && array_key_exists('localgrp',$rows[0]); if ($hasog)  $cols[] = array('origin (local)','s','ctr');
if ($IS_C) $cols[] = array('Cplus','s','ctr');
$cmarks = $IS_C ? cplus_load($CPLUS_JSON) : array();
$nrow = count($rows);
?><!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title><?= $esc($disp) ?> &middot; hDOM Database</title>
<style>
  :root {
    --bg: #fafaf7; --fg: #1d1d1b; --muted: #6b6b66; --grid: #e2e2dc; --card: #fff;
    --accent: #2f5fb3;
  }
  * { box-sizing: border-box; }
  body { margin: 0; padding: 32px 16px; background: var(--bg); color: var(--fg);
         font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
  .wrap { max-width: 1120px; margin: 0 auto; }
  a { color: var(--accent); text-decoration: none; }
  a:hover { text-decoration: underline; }
  .back { display: inline-block; margin-bottom: 14px; font-size: 13px; color: var(--muted); }
  .back:hover { color: var(--accent); text-decoration: none; }
  h1 { font-size: 24px; margin: 0 0 4px; display: flex; align-items: center; gap: 10px; }
  .dot { width: 12px; height: 12px; border-radius: 50%; flex-shrink: 0; }
  .sub { color: var(--muted); margin: 0; max-width: 760px; }
  .meta { display: flex; flex-wrap: wrap; align-items: baseline; gap: 6px 16px; margin: 14px 0 12px; color: var(--muted); font-size: 13px; }
  .meta b { font-size: 20px; color: var(--fg); font-variant-numeric: tabular-nums; }
  .meta .mark { color: #b8860b; font-weight: 600; }
  .tablewrap { overflow-x: auto; background: var(--card); border: 1px solid var(--grid); border-radius: 8px; }
  table { border-collapse: collapse; width: 100%; font-size: 13px; }
  th, td { padding: 7px 12px; border-bottom: 1px solid var(--grid); text-align: left; white-space: nowrap; }
  tbody tr:last-child td { border-bottom: 0; }
  th { font-weight: 600; color: var(--muted); font-size: 12px; cursor: pointer; user-select: none; }
  th:hover { color: var(--accent); }
  th[data-dir="asc"]::after { content: " \25B4"; }
  th[data-dir="desc"]::after { content: " \25BE"; }
  tbody tr:hover { background: #f4f6fb; }
  .idx { color: #9a9a94; }
  .num { text-align: right; font-variant-numeric: tabular-nums; }
  th.num { text-align: right; }
  .ctr { text-align: center; }
  td.mol a { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12.5px; }
</style>
</head>
<body>
<div class="wrap">
  <a class="back" href="../">&larr; <?= $esc($back) ?></a>
  <h1><?php if ($gcol !== ''): ?><span class="dot" style="background:<?= $esc($gcol) ?>"></span><?php endif; ?><?= $esc($disp) ?></h1>
<?php if ($gdesc !== ''): ?>  <p class="sub"><?= $esc($gdesc) ?></p>
<?php endif; ?>
  <div class="meta">
    <span><b><?= number_format($nrow) ?></b> <?= $nrow == 1 ? 'molecule' : 'molecules' ?></span>
    <a href="_summary.csv">_summary.csv</a>
<?php if ($IS_C): ?>    <span class="mark">&#9873; <span id="mvn"><?= count($cmarks) ?></span> marked &rarr; Cplus</span>
    <span>(open a molecule, press <b style="font-size:inherit">c</b> to mark; run move_cplus.py to move)</span>
<?php endif; ?>
  </div>
  <div class="tablewrap">
  <table id="t">
    <thead><tr><?php foreach ($cols as $k => $c): ?><th class="<?= $c[2] ?>" onclick="sortTable(<?= $k ?>,'<?= $c[1] ?>')"><?= $esc($c[0]) ?></th><?php endforeach; ?></tr></thead>
    <tbody>
<?php $i = 0; foreach ($rows as $r): $i++; $mol = $r['molecule'] ?? ''; if ($mol === '') continue; $v = 'index.php?file='.rawurlencode($mol);
      $hom = ($r['hom12'] ?? '') === 'Y'; $chg = (string)($r['change'] ?? ''); ?>
      <tr data-mol="<?= $esc($mol) ?>"><td class="idx"><?= $i ?></td><td class="mol"><a href="<?= $esc($v) ?>"><?= $esc($mol) ?></a></td><td class="num"><?= isset($LEN[$mol]) ? $esc($LEN[$mol]) : '' ?></td><td class="num"><?= $esc($r['score'] ?? '') ?></td><?php if (!$MALONLY): ?><td><?= $esc($r['local1'] ?? '') ?></td><?php endif; ?><td><?= $esc($r['mal1'] ?? '') ?></td><?php
      if ($hashom): ?><td class="ctr" style="font-weight:600;color:<?= $hom ? '#27ae60' : '#bbb' ?>"><?= $hom ? '&#10004;' : '&ndash;' ?></td><?php endif;
      if ($haschg): ?><td class="ctr" style="font-weight:600;color:<?= $chg === 'GAIN' ? '#27ae60' : ($chg === 'BREAK' ? '#e74c3c' : '#bbb') ?>"><?= $esc($chg) ?></td><?php endif;
      if ($hasog): ?><td class="ctr" style="font-weight:600"><?= $esc($r['localgrp'] ?? '') ?></td><?php endif;
      if ($IS_C): ?><td class="ctr cpm" data-mol="<?= $esc($mol) ?>" style="color:#b8860b;font-weight:600"><?= !empty($cmarks[$mol]) ? '&#9873;' : '' ?></td><?php endif; ?></tr>
<?php endforeach; ?>
    </tbody>
  </table>
  </div>
</div>
<script>
function sortTable(n,type){var t=document.getElementById('t'),tb=t.tBodies[0];var rows=Array.from(tb.rows);
var th=t.tHead.rows[0].cells[n];var dir=th.getAttribute('data-dir')==='asc'?'desc':'asc';
Array.from(t.tHead.rows[0].cells).forEach(function(c){if(c!==th)c.removeAttribute('data-dir');});th.setAttribute('data-dir',dir);
rows.sort(function(a,b){var x=a.cells[n].textContent.trim(),y=b.cells[n].textContent.trim();
if(type==='n'){x=parseFloat(x);y=parseFloat(y);if(isNaN(x))x=Infinity;if(isNaN(y))y=Infinity;}
if(x<y)return dir==='asc'?-1:1;if(x>y)return dir==='asc'?1:-1;return 0;});
rows.forEach(function(r){tb.appendChild(r);});renum();saveNavOrder();saveSort(n,type,dir);}
// The '#' column is a display position, not an identity: after any sort it reads 1..N in the
// order shown, so the number a row carries always matches where it sits on screen.
function renum(){var tb=document.getElementById('t').tBodies[0];
Array.from(tb.rows).forEach(function(r,i){if(r.cells[0])r.cells[0].textContent=i+1;});}
function navKey(){return 'navorder:'+location.pathname.replace(/index\.php$/,'').replace(/\/$/,'');}
function saveNavOrder(){try{var tb=document.getElementById('t').tBodies[0];
var order=Array.from(tb.rows).map(function(r){return r.getAttribute('data-mol');}).filter(Boolean);
localStorage.setItem(navKey(),JSON.stringify(order));}catch(e){}}saveNavOrder();
// The sort state (column, format, direction) is remembered per folder. Opening a molecule page
// and coming back with Up or the browser's back button redraws the list, and the sort last chosen
// is restored there.
function sortKey(){return 'sortstate:'+location.pathname.replace(/index\.php$/,'').replace(/\/$/,'');}
function saveSort(n,type,dir){try{localStorage.setItem(sortKey(),JSON.stringify({n:n,type:type,dir:dir}));}catch(e){}}
function restoreSort(){try{var s=JSON.parse(localStorage.getItem(sortKey())||'null');if(!s)return;
var t=document.getElementById('t');if(!t||!t.tHead)return;var th=t.tHead.rows[0].cells[s.n];if(!th)return;
// sortTable flips the current direction, so the opposite of the stored direction is seeded before calling it.
th.setAttribute('data-dir',s.dir==='asc'?'desc':'asc');sortTable(s.n,s.type);}catch(e){}}
restoreSort();
</script>
</body>
</html>
