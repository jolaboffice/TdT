<?php
// Header layout for the molecule view: the arrangement above the graphs that the hDOM_B1/B2
// viewers use. It is appended as CSS and JS before </body> instead of rewriting the generated
// view HTML. Only the header is taken; the profile and cc graphs and the h1 analysis table
// are left out.
function hdom_b2_layout() {
    return <<<'HTML'
<style>
/* Whitespace above the title. .container padding-top 30px plus the h1 top margin of 21.44px
   leaves 51px; 13px with no h1 margin is a quarter of that. The grey margin outside the
   container (body padding-top 40px) is unchanged. */
.container{padding-top:13px}
.title-row h1{margin-top:0}
/* Matching Results takes the width of the black data box above it. The per-column minima .mm
   sets from column 6 on are released: at the standard 1140px the columns are already wider than
   those minima, and in a narrower window the minima are what would hold the table past the box
   width. */
.matching-table.mm{width:100%;max-width:100%}
.mm td:nth-child(n+6),.mm th:nth-child(n+6){min-width:0}
.mm td:nth-child(12),.mm th:nth-child(12){min-width:0}
/* One set of arrows, on the title line, sized and rounded as in B2. */
.title-counter{display:inline-flex;align-items:center;gap:6px}
.title-counter .title-arrows{display:inline-flex;gap:5px;margin-left:8px}
.title-counter .nav-btn{width:30px;height:30px;font-size:16px;border-radius:5px;
  display:inline-flex;align-items:center;justify-content:center;padding:0}
</style>
<script>
(function(){
  /* The ref_map bar goes to the bottom of the page, where B2 keeps it. The id lets the toolbar
     attach after it. */
  var refbar=null;
  document.querySelectorAll('body > div').forEach(function(d){
    if(refbar) return;
    var s=d.querySelector(':scope > span');
    if(s && s.textContent.trim()==='ref_map:') refbar=d;
  });
  if(refbar){ refbar.id='ref-map-controls'; refbar.style.margin='20px 0 12px';
              document.body.appendChild(refbar); }

  /* The arrows move from the separate nav bar to the title counter and that bar is dropped,
     which also removes the duplicated n/n. */
  var nav=document.querySelector('.progress-bar-container + .nav-bar'),
      counter=document.querySelector('.title-row .title-counter');
  if(nav && counter){
    var arrows=document.createElement('span'); arrows.className='title-arrows';
    nav.querySelectorAll('.nav-btn').forEach(function(b){ arrows.appendChild(b); });
    counter.appendChild(arrows); nav.remove();
  }

  /* The Class/Mark/Green toolbar goes to the bottom, as in B2. It is moved rather than removed,
     so the controls keep working. */
  var bar=document.querySelector('.eval-toolbar');
  if(bar){
    var ctl=document.getElementById('ref-map-controls');
    if(ctl) ctl.after(bar); else document.body.appendChild(bar);
    bar.style.margin='12px 0 80px';
  }

  /* An h3.table-section-title is dropped only when it repeats the h2 immediately above it. */
  document.querySelectorAll('h3.table-section-title').forEach(function(h3){
    var p=h3.previousElementSibling;
    if(p && p.tagName==='H2' && p.textContent.trim()===h3.textContent.trim()) h3.remove();
  });
})();
</script>
HTML;
}
