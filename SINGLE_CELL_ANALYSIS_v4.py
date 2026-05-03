"""
=============================================================================
 SINGLE CELL NEIGHBOR ANALYSIS — 2G + 3G + 4G   v4
 PFE Project 2026 — 
=============================================================================
 CORRECTIONS v4 :
  ✅ FIX #1 — POSITIONS UNIFIÉES PAR SITE :
       Toutes les cellules 2G/3G/4G d'un même site sont placées AU MÊME POINT
       géographique (coordonnées du site). Chaque secteur est dessiné avec son
       azimuth réel depuis CE POINT commun (pas de décalage ox/oy).
       → Plus de points éparpillés par technologie pour le même site.

  ✅ FIX #2 — ZONES DE COUVERTURE CONCENTRIQUES VISIBLES :
       Les wedges (secteurs) sont tracés depuis le même point :
         • 2G = rayon 200 m  (couleur verte)
         • 3G = rayon 300 m  (couleur orange)
         • 4G = rayon 400 m  (couleur bleue)
       Comme chaque techno a un rayon différent, les 3 secteurs sont
       bien visibles en superposition.

  ✅ Fix TypeError LabelPlacement (compatible toutes versions QGIS)
  ✅ Cellule REF = étoile colorée bien visible
  ✅ Lignes voisins ROUGE / VERT / ORANGE avec HSR% sur chaque ligne
  ✅ Interface de sélection avec filtres
=============================================================================
"""

import os
import processing
import pandas as pd
from math import radians, degrees, sin, cos, asin, atan2, sqrt
from qgis.core import (
    QgsProject, QgsVectorLayer, QgsFeature,
    QgsGeometry, QgsPointXY, QgsField,
    QgsSymbol, QgsSimpleLineSymbolLayer,
    QgsSimpleMarkerSymbolLayer, QgsSimpleFillSymbolLayer,
    QgsCategorizedSymbolRenderer, QgsRendererCategory,
    QgsSingleSymbolRenderer, QgsProperty,
    QgsPalLayerSettings, QgsVectorLayerSimpleLabeling,
    QgsTextFormat, QgsTextBufferSettings
)
from PyQt5.QtCore    import QVariant, Qt
from PyQt5.QtGui     import QColor, QFont
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel,
    QComboBox, QPushButton, QGroupBox,
    QLineEdit
)

# ═══════════════════════════════════════════════════════════════════════════
#  PARAMÈTRES
# ═══════════════════════════════════════════════════════════════════════════
HSR_THRESH  = 97.0
DIST_4G     = 4000
DIST_3G     = 3000
DIST_2G     = 6000

# ── Rayons de couverture par technologie ─────────────────────────────────
RADIUS_4G   = 400   # mètres — Bleu
RADIUS_3G   = 300   # mètres — Orange
RADIUS_2G   = 200   # mètres — Vert

BEAM_WIDTH  = 120   # largeur du secteur (degrés)

# ── Couleurs par technologie ──────────────────────────────────────────────
COLOR_4G_FILL   = QColor(30,  144, 255, 90)   # Bleu Dodger, semi-transparent
COLOR_3G_FILL   = QColor(255, 140,   0, 90)   # Orange, semi-transparent
COLOR_2G_FILL   = QColor(50,  205,  50, 90)   # Vert lime, semi-transparent

COLOR_4G_POINT  = QColor(0,   100, 220)        # Bleu
COLOR_3G_POINT  = QColor(220, 100,   0)        # Orange foncé
COLOR_2G_POINT  = QColor(0,   160,  50)        # Vert

# ═══════════════════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════════════════
def hav(lon1, lat1, lon2, lat2):
    R = 6371000
    lo1,la1,lo2,la2 = map(radians,[lon1,lat1,lon2,lat2])
    a = sin((la2-la1)/2)**2 + cos(la1)*cos(la2)*sin((lo2-lo1)/2)**2
    return R * 2 * asin(sqrt(a))

def hsr_safe(att, fail):
    return 0.0 if att<=0 else round((att-fail)/att*100, 1)

def set_placement(lbl, name):
    """Compatible QGIS 3.16 → 3.36+ (enum renommé selon la version)"""
    for attr in [f"Placement.{name}", name, f"LabelPlacement.{name}"]:
        try:
            parts = attr.split(".")
            obj = QgsPalLayerSettings
            for p in parts:
                obj = getattr(obj, p)
            lbl.placement = obj
            return
        except AttributeError:
            continue
    lbl.placement = 0   # fallback OverPoint

def load_csv(path, label):
    if not os.path.exists(path):
        print(f"  ✗ MANQUANT: {os.path.basename(path)}")
        return None
    df = pd.read_csv(path, sep=";", encoding="utf-8-sig")
    df.columns = df.columns.str.strip()
    for c in df.select_dtypes("object").columns:
        df[c] = df[c].str.strip()
    print(f"  ✓ {label}: {len(df)} lignes")
    return df

def load_xlsx(path, label):
    if not os.path.exists(path):
        return None
    df = pd.read_excel(path)
    df.columns = df.columns.str.strip()
    for c in df.select_dtypes("object").columns:
        df[c] = df[c].str.strip()
    print(f"  ✓ {label}: {len(df)} lignes")
    return df

# ─────────────────────────────────────────────────────────────────────────
# FIX #1 : make_cell_map — on ne fait PLUS de décalage ox/oy
#   ox = lon  et  oy = lat  (position réelle du site)
#   Le secteur sera dessiné depuis ce point avec l'azimuth correct.
# ─────────────────────────────────────────────────────────────────────────
def make_cell_map(df, name_col, x_col, y_col, az_col, extras={}):
    """
    Construit un dictionnaire  cellid → {lon, lat, az, ox, oy, ...}
    ox et oy sont IDENTIQUES à lon/lat (pas de décalage).
    Cela garantit que toutes les technologies sont au même point géographique
    pour un même site, les secteurs étant différenciés par l'azimuth uniquement.
    """
    out = {}
    for _, r in df.iterrows():
        cid = str(r[name_col]).strip()
        try:
            lon = float(r[x_col])
            lat = float(r[y_col])
            az  = float(r[az_col]) if az_col and az_col in r.index else 0.0
        except:
            continue
        # ox/oy = position réelle du site (PAS de décalage)
        e = {"lon": lon, "lat": lat, "az": az, "ox": lon, "oy": lat}
        for k, col in extras.items():
            e[k] = r.get(col, "")
        out[cid] = e
    return out

def remove_layer(name):
    for lyr in QgsProject.instance().mapLayersByName(name):
        QgsProject.instance().removeMapLayer(lyr.id())

def determine_status(isn, hsr, att):
    if isn == 0: return "MISSING"
    if att > 0 and hsr >= HSR_THRESH: return "CONFIGURED"
    if att > 0 and hsr < HSR_THRESH:  return "POOR_HSR"
    return "MISSING"

def find_db(P, base):
    for sfx in [" (1).xlsx",".xlsx"]:
        fp = P+base+sfx
        if os.path.exists(fp): return fp
    return None

def make_label(field_expr, is_expr=False, size=7, placement="Line"):
    """Crée un QgsPalLayerSettings prêt à l'emploi."""
    txt = QgsTextFormat()
    txt.setFont(QFont("Arial", size, QFont.Bold))
    txt.setSize(size)
    buf = QgsTextBufferSettings()
    buf.setEnabled(True); buf.setSize(0.8); buf.setColor(QColor(255,255,255))
    txt.setBuffer(buf)
    lbl = QgsPalLayerSettings()
    lbl.setFormat(txt)
    lbl.fieldName    = field_expr
    lbl.isExpression = is_expr
    set_placement(lbl, placement)
    lbl.centroidWhole = True
    return lbl

# ═══════════════════════════════════════════════════════════════════════════
#  COUCHES QGIS
# ═══════════════════════════════════════════════════════════════════════════

def add_all_cell_points(name, cell_map, color, tech_label, ref_id=None):
    """
    Affiche TOUTES les cellules comme petits cercles (contexte réseau).
    La cellule REF est affichée en étoile plus grande avec son nom.
    Toutes les cellules sont placées à la position réelle du site (lon, lat).
    """
    remove_layer(name)
    vl = QgsVectorLayer("Point?crs=EPSG:4326", name, "memory")
    pr = vl.dataProvider()
    pr.addAttributes([
        QgsField("CELL_ID", QVariant.String),
        QgsField("IS_REF",  QVariant.Int),
        QgsField("AZIMUTH", QVariant.Double),
        QgsField("TECH",    QVariant.String),
    ])
    vl.updateFields()

    feats = []
    for cid, info in cell_map.items():
        try:
            f = QgsFeature()
            # FIX : utiliser lon/lat (position réelle), pas ox/oy
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(info["lon"], info["lat"])))
            is_ref = 1 if cid == ref_id else 0
            f.setAttributes([cid, is_ref, float(info.get("az", 0)), tech_label])
            feats.append(f)
        except: pass
    pr.addFeatures(feats)
    vl.updateExtents()

    # Renderer catégorisé : REF=étoile grande / autres=cercle petit
    cats = []
    sym_n = QgsSymbol.defaultSymbol(vl.geometryType())
    ml_n  = QgsSimpleMarkerSymbolLayer()
    ml_n.setShape(QgsSimpleMarkerSymbolLayer.Circle)
    ml_n.setColor(color); ml_n.setStrokeColor(color.darker(150))
    ml_n.setStrokeWidth(0.3); ml_n.setSize(3.0)
    sym_n.changeSymbolLayer(0, ml_n)
    cats.append(QgsRendererCategory(0, sym_n, f"Cellule {tech_label}"))

    sym_r = QgsSymbol.defaultSymbol(vl.geometryType())
    ml_r  = QgsSimpleMarkerSymbolLayer()
    ml_r.setShape(QgsSimpleMarkerSymbolLayer.Star)
    ml_r.setColor(QColor(255,220,0)); ml_r.setStrokeColor(QColor(180,120,0))
    ml_r.setStrokeWidth(0.8); ml_r.setSize(8.0)
    sym_r.changeSymbolLayer(0, ml_r)
    cats.append(QgsRendererCategory(1, sym_r, f"Cellule REF {tech_label} ★"))
    vl.setRenderer(QgsCategorizedSymbolRenderer("IS_REF", cats))

    # Label : nom complet sur la REF, rien sur les autres
    lbl = make_label(
        'if("IS_REF" = 1, "CELL_ID", \'\')',
        is_expr=True, size=9, placement="OverPoint"
    )
    vl.setLabeling(QgsVectorLayerSimpleLabeling(lbl))
    vl.setLabelsEnabled(True)

    QgsProject.instance().addMapLayer(vl)
    print(f"  → '{name}': {len(feats)} cellules  (REF=★ {ref_id})")
    return vl


def add_neighbor_lines(name, rows):
    """
    Lignes de la cellule REF vers chaque voisin.
    ROUGE  = MISSING    (IS_NEIGH=0)
    VERT   = CONFIGURED (IS_NEIGH=1, HSR ≥ 97%)
    ORANGE = POOR_HSR   (IS_NEIGH=1, HSR < 97%)
    Label = HSR% sur chaque ligne avec trafic.
    """
    remove_layer(name)
    if not rows:
        print(f"  (aucune ligne pour '{name}')")
        return None

    vl = QgsVectorLayer("LineString?crs=EPSG:4326", name, "memory")
    pr = vl.dataProvider()
    pr.addAttributes([
        QgsField("S_CELL",   QVariant.String),
        QgsField("T_CELL",   QVariant.String),
        QgsField("STATUS",   QVariant.String),
        QgsField("IS_NEIGH", QVariant.Int),
        QgsField("HO_ATT",   QVariant.Int),
        QgsField("HSR",      QVariant.Double),
        QgsField("DIST_M",   QVariant.Double),
        QgsField("RAT",      QVariant.String),
    ])
    vl.updateFields()

    feats = []
    for r in rows:
        try:
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(float(r["sox"]), float(r["soy"])),
                QgsPointXY(float(r["tox"]), float(r["toy"]))
            ]))
            f.setAttributes([r["sc"],r["tc"],r["status"],
                              int(r.get("isn",0)),int(r.get("att",0)),
                              float(r.get("hsr",0)),float(r.get("dist",0)),
                              r.get("rat","")])
            feats.append(f)
        except Exception as e:
            print(f"  Erreur feature: {e}")

    pr.addFeatures(feats)
    vl.updateExtents()

    color_cfg = {
        "MISSING":    (QColor(220,30,30),  2.8, "Voisin MANQUANT   (IS_NEIGH=0)"),
        "CONFIGURED": (QColor(30,180,30),  1.6, "Voisin CONFIGURÉ  (HSR ≥ 97%)"),
        "POOR_HSR":   (QColor(255,140,0),  2.2, "Voisin HSR FAIBLE (HSR < 97%)"),
    }
    cats = []
    for status,(col,w,lbl_txt) in color_cfg.items():
        sym = QgsSymbol.defaultSymbol(vl.geometryType())
        sl  = QgsSimpleLineSymbolLayer()
        sl.setColor(col); sl.setWidth(w)
        sym.changeSymbolLayer(0, sl)
        cats.append(QgsRendererCategory(status, sym, lbl_txt))
    vl.setRenderer(QgsCategorizedSymbolRenderer("STATUS", cats))

    lbl = make_label(
        'if("HSR" > 0, concat(to_string(round("HSR")), \'%\'), \'\')',
        is_expr=True, size=7, placement="Line"
    )
    vl.setLabeling(QgsVectorLayerSimpleLabeling(lbl))
    vl.setLabelsEnabled(True)

    QgsProject.instance().addMapLayer(vl)
    nm = sum(1 for r in rows if r["status"]=="MISSING")
    nc = sum(1 for r in rows if r["status"]=="CONFIGURED")
    np_ = sum(1 for r in rows if r["status"]=="POOR_HSR")
    print(f"  → '{name}': {len(feats)} lignes  "
          f"[🔴 {nm} manquants | 🟢 {nc} configurés | 🟠 {np_} HSR faible]")
    return vl


# ─────────────────────────────────────────────────────────────────────────
# FIX #2 : add_wedge_layer — rayons différents par technologie
#   Les wedges partent du point réel du site (lon, lat).
#   2G (200m vert) < 3G (300m orange) < 4G (400m bleu)
#   → Les 3 secteurs sont superposés et tous visibles car rayons croissants.
# ─────────────────────────────────────────────────────────────────────────
def add_wedge_layer(name, cell_map, fill_color, radius_m, ref_id=None):
    """
    Secteur d'antenne (wedge) pour TOUTES les cellules.
    - Le point de départ est la position réelle du site (lon, lat).
    - Le rayon varie selon la technologie : 2G=200m, 3G=300m, 4G=400m.
    - Les 3 secteurs (2G/3G/4G) seront superposés et visibles car
      leurs rayons sont différents + couleurs différentes.
    Si ref_id fourni : dessine uniquement le secteur de la cellule REF.
    """
    remove_layer(name)
    cells_to_draw = cell_map if ref_id is None else \
                    ({ref_id: cell_map[ref_id]} if ref_id in cell_map else {})
    if not cells_to_draw:
        print(f"  (pas de wedge pour '{name}')")
        return None

    radius_deg = radius_m / 111111.0

    tmp = QgsVectorLayer("Point?crs=EPSG:4326", "_tmp_wedge", "memory")
    tmp_pr = tmp.dataProvider()
    tmp_pr.addAttributes([
        QgsField("CELL_ID", QVariant.String),
        QgsField("AZIMUTH", QVariant.Double)
    ])
    tmp.updateFields()

    feats = []
    for cid, info in cells_to_draw.items():
        try:
            f = QgsFeature()
            # FIX : utiliser lon/lat (position réelle du site)
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(info["lon"], info["lat"])))
            f.setAttributes([cid, float(info.get("az", 0))])
            feats.append(f)
        except: pass

    tmp_pr.addFeatures(feats)
    tmp.updateExtents()

    params = {
        "INPUT":        tmp,
        "AZIMUTH":      QgsProperty.fromExpression('"AZIMUTH"'),
        "WIDTH":        BEAM_WIDTH / 2.0,
        "OUTER_RADIUS": radius_deg,
        "INNER_RADIUS": 0,
        "OUTPUT":       "memory:" + name,
    }
    res = processing.run("native:wedgebuffers", params)
    wl  = res["OUTPUT"]
    wl.setName(name)

    border = fill_color.darker(160)
    sym = QgsSymbol.defaultSymbol(wl.geometryType())
    fl  = QgsSimpleFillSymbolLayer()
    fl.setColor(fill_color)
    fl.setStrokeColor(border)
    fl.setStrokeWidth(0.8)
    sym.changeSymbolLayer(0, fl)
    wl.setRenderer(QgsSingleSymbolRenderer(sym))

    QgsProject.instance().addMapLayer(wl)
    print(f"  → Wedge '{name}': {len(feats)} secteurs  (rayon={radius_m}m)")
    return wl


# ═══════════════════════════════════════════════════════════════════════════
#  CHARGEMENT FICHIERS
# ═══════════════════════════════════════════════════════════════════════════
P = QgsProject.instance().readPath("./") + "/"
print("="*65)
print("  Chargement des fichiers...")
print(f"  Dossier: {P}")
print("="*65)

sho  = load_csv(P+"3G_SHO.csv",      "3G_SHO")
ifho = load_csv(P+"3G_IFHO.csv",     "3G_IFHO")
ho4g = load_csv(P+"4G_HO.csv",       "4G_HO")
ho2g = load_csv(P+"2G_HO.csv",       "2G_HO")
ref4g= load_csv(P+"4G_CELL_REF.csv", "4G_CELL_REF")
ref3g= load_csv(P+"3G_CELL_REF.csv", "3G_CELL_REF")
ref2g= load_csv(P+"2G_CELL_REF.csv", "2G_CELL_REF")

miss = [n for n,d in [("3G_SHO",sho),("3G_IFHO",ifho),("4G_HO",ho4g),
        ("2G_HO",ho2g),("4G_CELL_REF",ref4g),("3G_CELL_REF",ref3g),
        ("2G_CELL_REF",ref2g)] if d is None]
if miss:
    raise SystemExit(f"Fichiers manquants: {miss}")

_p4=find_db(P,"DB4G"); db4g=load_xlsx(_p4,"DB4G") if _p4 else None
_p3=find_db(P,"DB3G"); db3g=load_xlsx(_p3,"DB3G") if _p3 else None
_p2=find_db(P,"DB2G"); db2g=load_xlsx(_p2,"DB2G") if _p2 else None

# ── Cartes de cellules ────────────────────────────────────────────────────
# FIX : make_cell_map ne fait plus de décalage → ox=lon, oy=lat
m4g_ref = make_cell_map(ref4g,"CELLNAME","X","Y","AZIMUTH",{"site":"SITENAME"})
m3g_ref = make_cell_map(ref3g,"CELLNAME","X","Y","AZIMUTH",{"site":"SITENAME"})
m2g_ref = make_cell_map(ref2g,"CELLNAME","X","Y","AZIMUTH",{"site":"SITENAME"})

freq_col = "DL FREQ" if "DL FREQ" in ref3g.columns else "DL_FREQ"
if freq_col in ref3g.columns:
    for _,r in ref3g.iterrows():
        cid = str(r["CELLNAME"]).strip()
        if cid in m3g_ref:
            try: m3g_ref[cid]["freq"]=int(r[freq_col])
            except: pass

m4g_new={}
if db4g is not None:
    for _,r in db4g.iterrows():
        try:
            cid=str(r["EUtranCellFDDId"]).strip()
            lon,lat,az=float(r["LONG"]),float(r["LAT"]),float(r["AZIMUTH"])
            m4g_new[cid]={"lon":lon,"lat":lat,"az":az,"ox":lon,"oy":lat,
                           "site":str(r["Site_Name"]),"src":"DB"}
        except: pass

m3g_new={}
if db3g is not None:
    for _,r in db3g.iterrows():
        try:
            cid=str(r["UtranCellId"]).strip()
            lon,lat,az=float(r["LONG"]),float(r["LAT"]),float(r["AZIMUTH"])
            m3g_new[cid]={"lon":lon,"lat":lat,"az":az,"ox":lon,"oy":lat,
                           "site":str(r["3GSiteName"]),"freq":int(r["uarfcnDl"]),"src":"DB"}
        except: pass

m2g_new={}
if db2g is not None:
    sec_cnt={}
    for _,r in db2g.iterrows():
        try:
            site=str(r["SiteName"]).strip()
            sec=str(r.get("sector","")).strip()
            cid=f"{site}_{sec}" if sec else f"{site}_S{sec_cnt.get(site,0)+1}"
            sec_cnt[site]=sec_cnt.get(site,0)+1
            lon,lat,az=float(r["LONG"]),float(r["LAT"]),float(r["Azimuth"])
            m2g_new[cid]={"lon":lon,"lat":lat,"az":az,"ox":lon,"oy":lat,"site":site,"src":"DB"}
        except: pass

m4g_all={**m4g_ref,**m4g_new}
m3g_all={**m3g_ref,**m3g_new}
m2g_all={**m2g_ref,**m2g_new}

print(f"  4G: {len(m4g_all)} cellules  |  3G: {len(m3g_all)}  |  2G: {len(m2g_all)}")

# Listes pour les dropdowns
cells_4g = sorted(ho4g["S_CELL"].dropna().unique().tolist())
cells_3g = sorted(sho["S_CELL"].dropna().unique().tolist())
cells_2g = sorted(ho2g["S_CELL"].dropna().unique().tolist())


# ═══════════════════════════════════════════════════════════════════════════
#  FENÊTRE DE SÉLECTION
# ═══════════════════════════════════════════════════════════════════════════
class CellSelector(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Sélection des cellules de référence")
        self.setMinimumWidth(500)
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)

        lay = QVBoxLayout(self)
        lay.setSpacing(10)

        title = QLabel("📡  Choisissez une cellule de référence par technologie")
        title.setStyleSheet("font-size:13px;font-weight:bold;padding:6px;")
        title.setAlignment(Qt.AlignCenter)
        lay.addWidget(title)

        info = QLabel(
            f"  Seuil HSR : {HSR_THRESH}%   "
            f"[ 🟢 vert ≥ {HSR_THRESH}%  |  🟠 orange < {HSR_THRESH}%  |  🔴 rouge = manquant ]\n"
            f"  Rayons : 2G=200m (vert)  |  3G=300m (orange)  |  4G=400m (bleu)"
        )
        info.setStyleSheet("color:#444;font-size:11px;padding:2px 6px;")
        lay.addWidget(info)

        def make_group(title_txt, color, cells, attr_search, attr_combo):
            grp = QGroupBox(f"  {title_txt}")
            grp.setStyleSheet(f"QGroupBox{{font-weight:bold;color:{color};}}")
            v = QVBoxLayout(grp)
            search = QLineEdit()
            search.setPlaceholderText(f"Filtrer les cellules {title_txt[:2]}…")
            combo  = QComboBox()
            combo.setMaxVisibleItems(15)
            combo.addItems(cells)
            search.textChanged.connect(
                lambda t,c=combo,cl=cells: (c.clear(), c.addItems(
                    [x for x in cl if t.lower() in x.lower()])))
            v.addWidget(search); v.addWidget(combo)
            setattr(self, attr_search, search)
            setattr(self, attr_combo,  combo)
            return grp

        lay.addWidget(make_group("4G — Cellule de référence (Bleu  | 400m)",  "#0066cc", cells_4g, "s4","c4"))
        lay.addWidget(make_group("3G — Cellule de référence (Orange | 300m)", "#cc6600", cells_3g, "s3","c3"))
        lay.addWidget(make_group("2G — Cellule de référence (Vert   | 200m)", "#006633", cells_2g, "s2","c2"))

        btn_row = QHBoxLayout()
        btn_ok  = QPushButton("▶  Lancer l'analyse")
        btn_ok.setStyleSheet(
            "QPushButton{background:#0066cc;color:white;font-weight:bold;"
            "padding:9px 22px;border-radius:4px;font-size:13px;}"
            "QPushButton:hover{background:#0055aa;}")
        btn_ok.clicked.connect(self.accept)
        btn_cancel = QPushButton("Annuler")
        btn_cancel.setStyleSheet("padding:9px 14px;font-size:12px;")
        btn_cancel.clicked.connect(self.reject)
        btn_row.addStretch(); btn_row.addWidget(btn_ok); btn_row.addWidget(btn_cancel)
        lay.addLayout(btn_row)

    def cells(self):
        return self.c4.currentText(), self.c3.currentText(), self.c2.currentText()


dlg = CellSelector()
if not dlg.exec_():
    raise SystemExit("Annulé.")

REF_4G, REF_3G, REF_2G = dlg.cells()
print(f"\n  4G → {REF_4G}")
print(f"  3G → {REF_3G}")
print(f"  2G → {REF_2G}\n")


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSE 4G
# ═══════════════════════════════════════════════════════════════════════════
print("─"*65)
print(f"  4G : {REF_4G}")
g4_rows = []

if REF_4G in m4g_all:
    si = m4g_all[REF_4G]
    if "HO_FAIL" in ho4g.columns:
        ho4g["HO_HSR"] = ho4g.apply(lambda r: hsr_safe(r["HO_ATT"],r["HO_FAIL"]),axis=1)
    elif "HO_HSR" not in ho4g.columns:
        ho4g["HO_HSR"] = 0.0

    known = set()
    for _,r in ho4g[ho4g["S_CELL"]==REF_4G].iterrows():
        tc=r["T_CELL"]; isn=int(r["IS_NEIGH"]); att=int(r["HO_ATT"])
        hsr=float(r.get("HO_HSR",0)); dist=float(r["DIST"])
        ti=m4g_all.get(tc)
        if ti: tox,toy=ti["ox"],ti["oy"]
        else:
            try: tox,toy=float(r["T_X"]),float(r["T_Y"])
            except: continue
        g4_rows.append({"sc":REF_4G,"tc":tc,"status":determine_status(isn,hsr,att),
            "isn":isn,"att":att,"hsr":hsr,"dist":dist,"rat":"4G_HO",
            "sox":si["ox"],"soy":si["oy"],"tox":tox,"toy":toy})
        known.add(tc)

    for tc,ti in m4g_all.items():
        if tc==REF_4G or tc in known: continue
        d=hav(si["lon"],si["lat"],ti["lon"],ti["lat"])
        if d<=DIST_4G:
            g4_rows.append({"sc":REF_4G,"tc":tc,"status":"MISSING",
                "isn":0,"att":0,"hsr":0.0,"dist":round(d,1),"rat":"4G_SPATIAL",
                "sox":si["ox"],"soy":si["oy"],"tox":ti["ox"],"toy":ti["oy"]})
            known.add(tc)

    nm=sum(1 for r in g4_rows if r["status"]=="MISSING")
    nc=sum(1 for r in g4_rows if r["status"]=="CONFIGURED")
    np_=sum(1 for r in g4_rows if r["status"]=="POOR_HSR")
    print(f"  Total: {len(g4_rows)}  [🔴 {nm} | 🟢 {nc} | 🟠 {np_}]")
else:
    print(f"  ✗ '{REF_4G}' introuvable dans les données 4G")

# Ordre d'ajout : 2G en bas (rayon petit), 3G, puis 4G en haut (rayon grand)
# Pour la visibilité on ajoute d'abord le plus grand rayon (4G) en dessous
# puis 3G puis 2G au-dessus → ainsi les petits rayons ne sont pas masqués

# ─── COUCHES 4G ───────────────────────────────────────────────────────────
# Wedge 4G : bleu, rayon 400m — ajouté en premier (dessiné en dessous)
add_wedge_layer("Wedges_4G", m4g_all, COLOR_4G_FILL, RADIUS_4G, ref_id=None)
add_all_cell_points("Cell_Points_4G", m4g_all, COLOR_4G_POINT, "4G", ref_id=REF_4G)
add_neighbor_lines("REF_Neighbors_4G", g4_rows)


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSE 3G — SHO
# ═══════════════════════════════════════════════════════════════════════════
print("─"*65)
print(f"  3G SHO : {REF_3G}")
sho_rows=[]; ref_freq=0

if REF_3G in m3g_all:
    si3=m3g_all[REF_3G]; ref_freq=si3.get("freq",0)
    known=set()
    for _,r in sho[sho["S_CELL"]==REF_3G].iterrows():
        tc=r["T_CELL"]; isn=int(r["IS_NEIGH"]); att=int(r["HO_ATT"])
        hsr=float(r.get("SHO_HSR",0)); dist=float(r["DIST"])
        ti=m3g_all.get(tc)
        if ti: tox,toy=ti["ox"],ti["oy"]
        else:
            try: tox,toy=float(r["T_X"]),float(r["T_Y"])
            except: continue
        sho_rows.append({"sc":REF_3G,"tc":tc,"status":determine_status(isn,hsr,att),
            "isn":isn,"att":att,"hsr":hsr,"dist":dist,"rat":"3G_SHO",
            "sox":si3["ox"],"soy":si3["oy"],"tox":tox,"toy":toy})
        known.add(tc)
    for tc,ti in m3g_all.items():
        if tc==REF_3G or tc in known: continue
        if ref_freq and ti.get("freq",0) and ref_freq!=ti["freq"]: continue
        d=hav(si3["lon"],si3["lat"],ti["lon"],ti["lat"])
        if d<=DIST_3G:
            sho_rows.append({"sc":REF_3G,"tc":tc,"status":"MISSING",
                "isn":0,"att":0,"hsr":0.0,"dist":round(d,1),"rat":"3G_SHO_SPATIAL",
                "sox":si3["ox"],"soy":si3["oy"],"tox":ti["ox"],"toy":ti["oy"]})
            known.add(tc)
    print(f"  SHO: {len(sho_rows)}  [🔴{sum(1 for r in sho_rows if r['status']=='MISSING')} | 🟢{sum(1 for r in sho_rows if r['status']=='CONFIGURED')}]")
else:
    print(f"  ✗ '{REF_3G}' introuvable 3G")

# ── IFHO ─────────────────────────────────────────────────────────────────
print(f"  3G IFHO : {REF_3G}")
ifho_rows=[]
if REF_3G in m3g_all:
    known=set()
    for _,r in ifho[ifho["S_CELL"]==REF_3G].iterrows():
        tc=r["T_CELL"]
        sf=m3g_all.get(REF_3G,{}).get("freq",0)
        tf=m3g_all.get(tc,{}).get("freq",0)
        if sf and tf and sf==tf: continue
        isn=int(r["IS_NEIGH"]); att=int(r["HO_ATT"])
        hsr=float(r.get("IFHO_HSR",0)); dist=float(r["DIST"])
        ti=m3g_all.get(tc)
        if ti: tox,toy=ti["ox"],ti["oy"]
        else:
            try: tox,toy=float(r["T_X"]),float(r["T_Y"])
            except: continue
        ifho_rows.append({"sc":REF_3G,"tc":tc,"status":determine_status(isn,hsr,att),
            "isn":isn,"att":att,"hsr":hsr,"dist":dist,"rat":"3G_IFHO",
            "sox":si3["ox"],"soy":si3["oy"],"tox":tox,"toy":toy})
        known.add(tc)
    for tc,ti in m3g_all.items():
        if tc==REF_3G or tc in known: continue
        if ref_freq and ti.get("freq",0) and ref_freq==ti["freq"]: continue
        d=hav(si3["lon"],si3["lat"],ti["lon"],ti["lat"])
        if d<=DIST_3G:
            ifho_rows.append({"sc":REF_3G,"tc":tc,"status":"MISSING",
                "isn":0,"att":0,"hsr":0.0,"dist":round(d,1),"rat":"3G_IFHO_SPATIAL",
                "sox":si3["ox"],"soy":si3["oy"],"tox":ti["ox"],"toy":ti["oy"]})
    print(f"  IFHO: {len(ifho_rows)}  [🔴{sum(1 for r in ifho_rows if r['status']=='MISSING')} | 🟢{sum(1 for r in ifho_rows if r['status']=='CONFIGURED')}]")

# ─── COUCHES 3G ───────────────────────────────────────────────────────────
# Wedge 3G : orange, rayon 300m
add_wedge_layer("Wedges_3G", m3g_all, COLOR_3G_FILL, RADIUS_3G, ref_id=None)
add_all_cell_points("Cell_Points_3G", m3g_all, COLOR_3G_POINT, "3G", ref_id=REF_3G)
add_neighbor_lines("REF_Neighbors_3G_SHO",  sho_rows)
add_neighbor_lines("REF_Neighbors_3G_IFHO", ifho_rows)


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSE 2G
# ═══════════════════════════════════════════════════════════════════════════
print("─"*65)
print(f"  2G : {REF_2G}")
g2_rows=[]

if REF_2G in m2g_all:
    si2=m2g_all[REF_2G]
    ho2g["OUT_HSR"]=ho2g.apply(lambda r: hsr_safe(r["OUT_ATT"],r["OUT_FAIL"]),axis=1)
    ho2g["IN_HSR"] =ho2g.apply(lambda r: hsr_safe(r["IN_ATT"], r["IN_FAIL"]), axis=1)
    known=set()
    for _,r in ho2g[ho2g["S_CELL"]==REF_2G].iterrows():
        tc=r["T_CELL"]; isn=int(r["IS_NEIGH"])
        oa,oh=int(r["OUT_ATT"]),float(r["OUT_HSR"])
        ia,ih=int(r["IN_ATT"]), float(r["IN_HSR"])
        dist=float(r["DIST"])
        hsr=oh if oa>0 else ih; att=oa if oa>0 else ia
        ti=m2g_all.get(tc)
        if ti: tox,toy=ti["ox"],ti["oy"]
        else:
            try: tox,toy=float(r["T_X"]),float(r["T_Y"])
            except: continue
        g2_rows.append({"sc":REF_2G,"tc":tc,"status":determine_status(isn,hsr,att),
            "isn":isn,"att":att,"hsr":hsr,"dist":dist,"rat":"2G_HO",
            "sox":si2["ox"],"soy":si2["oy"],"tox":tox,"toy":toy})
        known.add(tc)
    for tc,ti in m2g_all.items():
        if tc==REF_2G or tc in known: continue
        d=hav(si2["lon"],si2["lat"],ti["lon"],ti["lat"])
        if d<=DIST_2G:
            g2_rows.append({"sc":REF_2G,"tc":tc,"status":"MISSING",
                "isn":0,"att":0,"hsr":0.0,"dist":round(d,1),"rat":"2G_SPATIAL",
                "sox":si2["ox"],"soy":si2["oy"],"tox":ti["ox"],"toy":ti["oy"]})
            known.add(tc)
    nm=sum(1 for r in g2_rows if r["status"]=="MISSING")
    nc=sum(1 for r in g2_rows if r["status"]=="CONFIGURED")
    np_=sum(1 for r in g2_rows if r["status"]=="POOR_HSR")
    print(f"  Total: {len(g2_rows)}  [🔴 {nm} | 🟢 {nc} | 🟠 {np_}]")
else:
    print(f"  ✗ '{REF_2G}' introuvable 2G")

# ─── COUCHES 2G ───────────────────────────────────────────────────────────
# Wedge 2G : vert, rayon 200m — ajouté en dernier (dessiné au-dessus)
add_wedge_layer("Wedges_2G", m2g_all, COLOR_2G_FILL, RADIUS_2G, ref_id=None)
add_all_cell_points("Cell_Points_2G", m2g_all, COLOR_2G_POINT, "2G", ref_id=REF_2G)
add_neighbor_lines("REF_Neighbors_2G", g2_rows)


# ── Rafraîchissement carte ────────────────────────────────────────────────
from qgis.utils import iface
iface.mapCanvas().refresh()

# ═══════════════════════════════════════════════════════════════════════════
print()
print("="*65)
print("  ✅  ANALYSE TERMINÉE")
print("="*65)
print(f"  REF 4G : {REF_4G}  →  {len(g4_rows)} voisins")
print(f"  REF 3G : {REF_3G}  →  SHO:{len(sho_rows)}  IFHO:{len(ifho_rows)}")
print(f"  REF 2G : {REF_2G}  →  {len(g2_rows)} voisins")
print()
print("  COUCHES CRÉÉES (ordre carte) :")
print("  ├─ Wedges_4G  → secteurs bleu   (rayon 400m)")
print("  ├─ Wedges_3G  → secteurs orange (rayon 300m)")
print("  ├─ Wedges_2G  → secteurs vert   (rayon 200m)  ← au-dessus")
print("  ├─ Cell_Points_4G / 3G / 2G → ● toutes cellules  ★ REF (étoile jaune)")
print("  ├─ REF_Neighbors_4G / 2G   → lignes voisins")
print("  ├─ REF_Neighbors_3G_SHO    → lignes voisins SHO")
print("  └─ REF_Neighbors_3G_IFHO   → lignes voisins IFHO")
print()
print("  LÉGENDE WEDGES (zones de couverture) :")
print("  🔵 BLEU   (400m) = secteur 4G")
print("  🟠 ORANGE (300m) = secteur 3G")
print("  🟢 VERT   (200m) = secteur 2G")
print()
print("  LÉGENDE LIGNES :")
print("  🔴 ROUGE  = voisin MANQUANT  (IS_NEIGH=0)")
print("  🟢 VERT   = voisin CONFIGURÉ (HSR ≥ 97%)")
print("  🟠 ORANGE = voisin en place  (HSR < 97%)")
print("  Label sur chaque ligne = HSR%")
print("="*65)
print()
print("  💡 ASTUCE : Dans le panneau des couches, placez les Wedges")
print("     dans cet ordre (haut→bas) pour une meilleure visibilité :")
print("     Wedges_2G → Wedges_3G → Wedges_4G")
print("="*65)
