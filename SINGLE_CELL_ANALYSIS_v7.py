"""
=============================================================================
 SINGLE CELL NEIGHBOR ANALYSIS — 2G + 3G + 4G   v4.1  (PRO)
 PFE Project 2026 — Tunisie Telecom
=============================================================================
 NOUVEAUTÉS v4 :
  ✅ Zoom automatique sur la cellule de référence après analyse
  ✅ Barre de progression pendant l'exécution
  ✅ Export automatique rapport .txt dans le dossier projet
  ✅ Épaisseur des lignes proportionnelle au trafic HO_ATT
  ✅ Vérification des colonnes avant traitement (plus de crash)
  ✅ Log complet écrit dans analysis_log.txt
  ✅ Bloc CONFIG unique en haut du script
  ✅ Popup résumé final avec statistiques complètes
 NOUVEAUTÉS v4.1 :
  ✅ Exactement 3 wedges par technologie (3 secteurs du même site REF)
  ✅ Même rayon pour tous les secteurs (RADIUS_UNIQUE)
  ✅ Secteur REF mis en évidence (contour plus épais)
=============================================================================
"""

import os, sys, logging
import processing
import pandas as pd
from datetime import datetime
from math import radians, degrees, sin, cos, asin, atan2, sqrt
from qgis.core import (
    QgsProject, QgsVectorLayer, QgsFeature,
    QgsGeometry, QgsPointXY, QgsField,
    QgsSymbol, QgsSimpleLineSymbolLayer,
    QgsSimpleMarkerSymbolLayer, QgsSimpleFillSymbolLayer,
    QgsCategorizedSymbolRenderer, QgsRendererCategory,
    QgsSingleSymbolRenderer, QgsProperty,
    QgsPalLayerSettings, QgsVectorLayerSimpleLabeling,
    QgsTextFormat, QgsTextBufferSettings,
    QgsRectangle, QgsCoordinateReferenceSystem
)
from PyQt5.QtCore    import QVariant, Qt
from PyQt5.QtGui     import QColor, QFont
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel,
    QComboBox, QPushButton, QGroupBox,
    QLineEdit, QProgressDialog, QMessageBox,
    QApplication, QSizePolicy
)

# ═══════════════════════════════════════════════════════════════════════════
#  ██████  BLOC CONFIG — MODIFIEZ ICI UNIQUEMENT  ██████
# ═══════════════════════════════════════════════════════════════════════════

CONFIG = {
    # ── Seuil HSR : voisin CONFIGURÉ si HSR ≥ ce seuil (%) ──────────────
    "HSR_THRESH"  : 97.0,

    # ── Distance max recherche voisins spatiaux (mètres) ─────────────────
    "DIST_4G"     : 4000,
    "DIST_3G"     : 3000,
    "DIST_2G"     : 6000,

    # ── Décalage visuel du point cellule (mètres) ────────────────────────
    "OFFSET_M"    : 120,

    # ── Rayons wedge distincts par technologie (mètres) ──────────────────
    "RADIUS_4G"   : 150,   # 4G = 150 m
    "RADIUS_3G"   : 100,   # 3G = 100 m
    "RADIUS_2G"   :  50,   # 2G =  50 m

    # ── Largeur du secteur antenne (degrés) ──────────────────────────────
    "BEAM_WIDTH"  : 120,

    # ── Épaisseur min/max des lignes voisins (selon trafic HO_ATT) ───────
    "LINE_W_MIN"  : 0.8,
    "LINE_W_MAX"  : 4.0,

    # ── Noms des fichiers CSV (dans le dossier du projet .qgz) ───────────
    "FILE_3G_SHO"      : "3G_SHO.csv",
    "FILE_3G_IFHO"     : "3G_IFHO.csv",
    "FILE_4G_HO"       : "4G_HO.csv",
    "FILE_2G_HO"       : "2G_HO.csv",
    "FILE_4G_CELL_REF" : "4G_CELL_REF.csv",
    "FILE_3G_CELL_REF" : "3G_CELL_REF.csv",
    "FILE_2G_CELL_REF" : "2G_CELL_REF.csv",
    "FILE_DB4G"        : "DB4G",
    "FILE_DB3G"        : "DB3G",
    "FILE_DB2G"        : "DB2G",
}

# ═══════════════════════════════════════════════════════════════════════════
#  LOGGER — console + fichier
# ═══════════════════════════════════════════════════════════════════════════
P = QgsProject.instance().readPath("./") + "/"
LOG_FILE = P + f"analysis_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

logging.basicConfig(
    level    = logging.INFO,
    format   = "%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt  = "%H:%M:%S",
    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, encoding="utf-8"),
    ]
)
log = logging.getLogger("CellAnalysis")
log.info("="*65)
log.info("  SINGLE CELL NEIGHBOR ANALYSIS  v4")
log.info(f"  Dossier projet : {P}")
log.info(f"  Log écrit dans : {LOG_FILE}")
log.info("="*65)

# Alias CONFIG
HSR_THRESH = CONFIG["HSR_THRESH"]
DIST_4G    = CONFIG["DIST_4G"]
DIST_3G    = CONFIG["DIST_3G"]
DIST_2G    = CONFIG["DIST_2G"]
OFFSET_M   = CONFIG["OFFSET_M"]
RADIUS_4G  = CONFIG["RADIUS_4G"]   # 150 m
RADIUS_3G  = CONFIG["RADIUS_3G"]   # 100 m
RADIUS_2G  = CONFIG["RADIUS_2G"]   #  50 m
BEAM_WIDTH = CONFIG["BEAM_WIDTH"]
LINE_W_MIN = CONFIG["LINE_W_MIN"]
LINE_W_MAX = CONFIG["LINE_W_MAX"]

# ═══════════════════════════════════════════════════════════════════════════
#  FONCTIONS UTILITAIRES
# ═══════════════════════════════════════════════════════════════════════════

def hav(lon1, lat1, lon2, lat2):
    """Distance haversine en mètres."""
    R = 6371000
    lo1,la1,lo2,la2 = map(radians,[lon1,lat1,lon2,lat2])
    a = sin((la2-la1)/2)**2 + cos(la1)*cos(la2)*sin((lo2-lo1)/2)**2
    return R * 2 * asin(sqrt(a))

def offset_cell(lon, lat, az_deg, dist_m=None):
    """Décale un point de dist_m mètres dans la direction de l'azimut."""
    if dist_m is None: dist_m = OFFSET_M
    R = 6371000
    az = radians(az_deg)
    la1,lo1 = radians(lat),radians(lon)
    d_r = dist_m / R
    la2 = asin(sin(la1)*cos(d_r)+cos(la1)*sin(d_r)*cos(az))
    lo2 = lo1 + atan2(sin(az)*sin(d_r)*cos(la1), cos(d_r)-sin(la1)*sin(la2))
    return degrees(lo2), degrees(la2)

def hsr_safe(att, fail):
    """Calcule le HSR% sans division par zéro."""
    return 0.0 if att <= 0 else round((att - fail) / att * 100, 1)

def check_columns(df, required_cols, file_label):
    """Vérifie que toutes les colonnes nécessaires existent."""
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        log.error(f"  [{file_label}] Colonnes MANQUANTES : {missing}")
        log.info (f"  Colonnes disponibles : {df.columns.tolist()}")
        raise SystemExit(f"Colonnes manquantes dans {file_label}: {missing}")
    return True

def load_csv(path, label, required_cols=None):
    """Charge un CSV avec vérification de colonnes."""
    if not os.path.exists(path):
        log.error(f"  FICHIER MANQUANT : {os.path.basename(path)}")
        return None
    df = pd.read_csv(path, sep=";", encoding="utf-8-sig")
    df.columns = df.columns.str.strip()
    for c in df.select_dtypes("object").columns:
        df[c] = df[c].str.strip()
    if required_cols:
        check_columns(df, required_cols, label)
    log.info(f"  ✓ {label}: {len(df)} lignes  |  colonnes: {df.columns.tolist()}")
    return df

def load_xlsx(path, label):
    """Charge un fichier Excel."""
    if not path or not os.path.exists(path):
        return None
    df = pd.read_excel(path)
    df.columns = df.columns.str.strip()
    for c in df.select_dtypes("object").columns:
        df[c] = df[c].str.strip()
    log.info(f"  ✓ {label}: {len(df)} lignes")
    return df

def find_xlsx(P, base):
    """Cherche un .xlsx avec ou sans suffixe (1)."""
    for sfx in [" (1).xlsx", ".xlsx"]:
        fp = P + base + sfx
        if os.path.exists(fp):
            return fp
    return None

def make_cell_map(df, name_col, x_col, y_col, az_col, extras={}):
    """Construit un dict {cell_id: {lon, lat, az, ox, oy, ...}}"""
    out = {}
    for _, r in df.iterrows():
        cid = str(r[name_col]).strip()
        try:
            lon=float(r[x_col]); lat=float(r[y_col])
            az=float(r[az_col]) if az_col and az_col in r.index else 0.0
        except:
            continue
        ox,oy = offset_cell(lon,lat,az)
        e = {"lon":lon,"lat":lat,"az":az,"ox":ox,"oy":oy}
        for k,col in extras.items():
            e[k] = r.get(col,"")
        out[cid] = e
    return out

def remove_layer(name):
    """Supprime une couche existante du même nom."""
    for lyr in QgsProject.instance().mapLayersByName(name):
        QgsProject.instance().removeMapLayer(lyr.id())

def determine_status(isn, hsr, att):
    """Détermine le statut d'un voisin."""
    if isn == 0:                          return "MISSING"
    if att > 0 and hsr >= HSR_THRESH:    return "CONFIGURED"
    if att > 0 and hsr <  HSR_THRESH:    return "POOR_HSR"
    return "MISSING"

def set_label_placement(lbl, mode="Line"):
    """Compatible QGIS 3.16 → 3.36+ pour le placement des labels."""
    try:
        lbl.placement = getattr(QgsPalLayerSettings.Placement, mode)
    except AttributeError:
        try:
            lbl.placement = getattr(QgsPalLayerSettings, mode)
        except AttributeError:
            lbl.placement = 2 if mode == "Line" else 0

def make_label(field_expr, is_expr=False, size=7, placement="Line",
               color=QColor(0,0,0)):
    """Crée un QgsPalLayerSettings configuré."""
    txt = QgsTextFormat()
    txt.setFont(QFont("Arial", size, QFont.Bold))
    txt.setSize(size)
    txt.setColor(color)
    buf = QgsTextBufferSettings()
    buf.setEnabled(True); buf.setSize(0.9); buf.setColor(QColor(255,255,255))
    txt.setBuffer(buf)
    lbl = QgsPalLayerSettings()
    lbl.setFormat(txt)
    lbl.fieldName    = field_expr
    lbl.isExpression = is_expr
    set_label_placement(lbl, placement)
    lbl.centroidWhole = True
    return lbl

def att_to_width(att, att_max):
    """Convertit HO_ATT en épaisseur de ligne (proportionnel)."""
    if att_max <= 0 or att <= 0:
        return LINE_W_MIN
    ratio = min(att / att_max, 1.0)
    return round(LINE_W_MIN + ratio * (LINE_W_MAX - LINE_W_MIN), 2)

def zoom_to_cell(cell_map, ref_id, margin_deg=0.02):
    """Zoom automatique sur la cellule de référence."""
    if ref_id not in cell_map:
        return
    info = cell_map[ref_id]
    lon, lat = info["lon"], info["lat"]
    rect = QgsRectangle(lon-margin_deg, lat-margin_deg,
                        lon+margin_deg, lat+margin_deg)
    from qgis.utils import iface
    iface.mapCanvas().setExtent(rect)
    iface.mapCanvas().refresh()

# ═══════════════════════════════════════════════════════════════════════════
#  COUCHES QGIS
# ═══════════════════════════════════════════════════════════════════════════

def add_all_cell_points(name, cell_map, color, tech_label, ref_id=None,
                        freq_map=None):
    """
    Toutes les cellules en cercles colorés + label  'NOM_CELLULE  [BANDE]'
    sur TOUTES les cellules.
      • Cellule normale  : cercle coloré, label taille 7 (noir, buffer blanc)
      • Cellule REF      : étoile jaune, label taille 9 (couleur tech, buffer blanc)
    freq_map : dict {cell_id: "L800"/"U900"/"G900"/...}  (peut être None)
    """
    if freq_map is None:
        freq_map = {}
    remove_layer(name)
    vl = QgsVectorLayer("Point?crs=EPSG:4326", name, "memory")
    pr = vl.dataProvider()
    pr.addAttributes([
        QgsField("CELL_ID",   QVariant.String),
        QgsField("IS_REF",    QVariant.Int),
        QgsField("AZIMUTH",   QVariant.Double),
        QgsField("TECH",      QVariant.String),
        QgsField("FREQ_BAND", QVariant.String),   # ← nouveau champ bande
        QgsField("LBL",       QVariant.String),   # ← "NOM  [BANDE]" pré-calculé
    ])
    vl.updateFields()

    feats = []
    for cid, info in cell_map.items():
        try:
            band = freq_map.get(cid.strip(), "?")
            lbl  = f"{cid}  [{band}]"
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPointXY(
                QgsPointXY(info["ox"], info["oy"])))
            f.setAttributes([
                cid,
                1 if cid == ref_id else 0,
                float(info.get("az", 0)),
                tech_label,
                band,
                lbl,
            ])
            feats.append(f)
        except:
            pass
    pr.addFeatures(feats)
    vl.updateExtents()

    # ── Renderer : cercle normal / étoile REF ────────────────────────────
    cats = []
    sym_n = QgsSymbol.defaultSymbol(vl.geometryType())
    ml_n  = QgsSimpleMarkerSymbolLayer()
    ml_n.setShape(QgsSimpleMarkerSymbolLayer.Circle)
    ml_n.setColor(color)
    ml_n.setStrokeColor(color.darker(150))
    ml_n.setStrokeWidth(0.3)
    ml_n.setSize(3.0)
    sym_n.changeSymbolLayer(0, ml_n)
    cats.append(QgsRendererCategory(0, sym_n, f"Cellule {tech_label}"))

    sym_r = QgsSymbol.defaultSymbol(vl.geometryType())
    ml_r  = QgsSimpleMarkerSymbolLayer()
    ml_r.setShape(QgsSimpleMarkerSymbolLayer.Star)
    ml_r.setColor(QColor(255, 220, 0))
    ml_r.setStrokeColor(QColor(160, 100, 0))
    ml_r.setStrokeWidth(0.9)
    ml_r.setSize(9.0)
    sym_r.changeSymbolLayer(0, ml_r)
    cats.append(QgsRendererCategory(1, sym_r, f"★ Cellule REF {tech_label}"))
    vl.setRenderer(QgsCategorizedSymbolRenderer("IS_REF", cats))

    # ── Labels : toutes cellules → "NOM  [BANDE]" ────────────────────────
    # REF : taille 9, couleur technologie  |  Autres : taille 7, noir
    # On utilise l'expression data-defined sur la taille du texte
    txt = QgsTextFormat()
    txt.setFont(QFont("Arial", 7, QFont.Bold))
    txt.setSize(7)
    txt.setColor(QColor(30, 30, 30))
    # Taille data-defined : 9 pour REF, 7 pour les autres
    txt.dataDefinedProperties().setProperty(
        QgsTextFormat.Property.Size
        if hasattr(QgsTextFormat, "Property") else 4,   # 4 = Size property index
        QgsProperty.fromExpression('if("IS_REF" = 1, 9, 7)')
    )
    # Couleur data-defined : couleur tech pour REF, noir pour les autres
    _rgb = f"{color.red()},{color.green()},{color.blue()}"
    txt.dataDefinedProperties().setProperty(
        QgsTextFormat.Property.Color
        if hasattr(QgsTextFormat, "Property") else 1,   # 1 = Color property index
        QgsProperty.fromExpression(
            f'if("IS_REF" = 1, color_rgb({_rgb}), color_rgb(30,30,30))'
        )
    )
    buf = QgsTextBufferSettings()
    buf.setEnabled(True)
    buf.setSize(1.0)
    buf.setColor(QColor(255, 255, 255))
    txt.setBuffer(buf)

    lbl = QgsPalLayerSettings()
    lbl.setFormat(txt)
    lbl.fieldName    = '"LBL"'          # champ pré-calculé "NOM  [BANDE]"
    lbl.isExpression = True
    set_label_placement(lbl, "OverPoint")
    lbl.centroidWhole = True

    vl.setLabeling(QgsVectorLayerSimpleLabeling(lbl))
    vl.setLabelsEnabled(True)

    QgsProject.instance().addMapLayer(vl)
    log.info(f"  → '{name}': {len(feats)} cellules  (REF=★ {ref_id})")
    return vl


def add_neighbor_lines(name, rows):
    """
    Lignes vers les voisins.
    ROUGE  = MISSING    | VERT = CONFIGURED | ORANGE = POOR_HSR
    Épaisseur ∝ trafic HO_ATT.
    Label HSR% coloré (vert/rouge) centré sur la ligne.
    """
    remove_layer(name)
    if not rows:
        log.info(f"  (aucune ligne pour '{name}')")
        return None

    # Calcul épaisseur proportionnelle au trafic
    att_values = [r.get("att",0) for r in rows if r["status"] != "MISSING"]
    att_max    = max(att_values) if att_values else 1

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
        QgsField("LINE_W",   QVariant.Double),   # épaisseur calculée
    ])
    vl.updateFields()

    feats = []
    for r in rows:
        try:
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(float(r["sox"]),float(r["soy"])),
                QgsPointXY(float(r["tox"]),float(r["toy"]))
            ]))
            w = att_to_width(r.get("att",0), att_max) if r["status"]!="MISSING" else LINE_W_MIN
            f.setAttributes([
                r["sc"],r["tc"],r["status"],
                int(r.get("isn",0)), int(r.get("att",0)),
                float(r.get("hsr",0)), float(r.get("dist",0)),
                r.get("rat",""), float(w)
            ])
            feats.append(f)
        except Exception as e:
            log.warning(f"  Feature ignorée : {e}")

    pr.addFeatures(feats)
    vl.updateExtents()

    # Renderer 3 catégories
    color_cfg = {
        "MISSING":    (QColor(220, 30, 30), "Voisin MANQUANT   (IS_NEIGH=0)"),
        "CONFIGURED": (QColor( 30,180, 30), "Voisin CONFIGURÉ  (HSR ≥ 97%)"),
        "POOR_HSR":   (QColor(255,140,  0), "Voisin HSR FAIBLE (HSR < 97%)"),
    }
    cats = []
    for status,(col,lbl_txt) in color_cfg.items():
        sym = QgsSymbol.defaultSymbol(vl.geometryType())
        sl  = QgsSimpleLineSymbolLayer()
        sl.setColor(col)
        # Épaisseur data-defined via expression sur le champ LINE_W
        sl.setDataDefinedProperty(
            QgsSimpleLineSymbolLayer.PropertyStrokeWidth,
            QgsProperty.fromExpression(
                f'if("STATUS" = \'MISSING\', {LINE_W_MIN}, "LINE_W")'
            )
        )
        sym.changeSymbolLayer(0, sl)
        cats.append(QgsRendererCategory(status, sym, lbl_txt))
    vl.setRenderer(QgsCategorizedSymbolRenderer("STATUS", cats))

    # Label HSR% coloré (expression couleur selon seuil)
    lbl = make_label(
        'if("HSR" > 0, concat(to_string(round("HSR")), \'%\'), \'\')',
        is_expr=True, size=7, placement="Line"
    )
    vl.setLabeling(QgsVectorLayerSimpleLabeling(lbl))
    vl.setLabelsEnabled(True)

    QgsProject.instance().addMapLayer(vl)
    nm = sum(1 for r in rows if r["status"]=="MISSING")
    nc = sum(1 for r in rows if r["status"]=="CONFIGURED")
    np_= sum(1 for r in rows if r["status"]=="POOR_HSR")
    log.info(f"  → '{name}': {len(feats)} lignes  "
             f"[🔴 {nm} manquants | 🟢 {nc} configurés | 🟠 {np_} HSR faible]")
    return vl


def add_wedge_layer(name, ref_id, cell_map, fill_color, radius_m,
                    beam_width=None):
    """
    Crée un wedge pour TOUTES les cellules de la carte (tous les sites).
    Même rayon pour tous (radius_m).
    Le secteur REF est mis en évidence avec un contour plus épais.
    """
    if beam_width is None:
        beam_width = BEAM_WIDTH
    remove_layer(name)

    if not cell_map:
        log.warning(f"  Aucune cellule pour '{name}'")
        return None

    radius_deg = radius_m / 111111.0

    tmp    = QgsVectorLayer("Point?crs=EPSG:4326", "_tmp_wedge", "memory")
    tmp_pr = tmp.dataProvider()
    tmp_pr.addAttributes([
        QgsField("CELL_ID", QVariant.String),
        QgsField("AZIMUTH", QVariant.Double),
        QgsField("IS_REF",  QVariant.Int),
    ])
    tmp.updateFields()

    feats = []
    for cid, info in cell_map.items():
        try:
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPointXY(
                QgsPointXY(info["lon"], info["lat"])))
            f.setAttributes([cid, float(info.get("az", 0)),
                              1 if cid == ref_id else 0])
            feats.append(f)
        except Exception as e:
            log.warning(f"  Wedge ignorée ({cid}): {e}")

    tmp_pr.addFeatures(feats)
    tmp.updateExtents()

    params = {
        "INPUT":        tmp,
        "AZIMUTH":      QgsProperty.fromExpression('"AZIMUTH"'),
        "WIDTH":        beam_width / 2.0,
        "OUTER_RADIUS": radius_deg,
        "INNER_RADIUS": 0,
        "OUTPUT":       "memory:" + name,
    }
    res = processing.run("native:wedgebuffers", params)
    wl  = res["OUTPUT"]
    wl.setName(name)

    border = fill_color.darker(170)

    # Secteurs normaux (tous les sites)
    sym_n = QgsSymbol.defaultSymbol(wl.geometryType())
    fl_n  = QgsSimpleFillSymbolLayer()
    fl_n.setColor(fill_color)
    fl_n.setStrokeColor(border)
    fl_n.setStrokeWidth(0.4)
    sym_n.changeSymbolLayer(0, fl_n)

    # Secteur REF — contour épais + plus opaque pour se démarquer
    ref_fill = QColor(fill_color)
    ref_fill.setAlpha(min(fill_color.alpha() + 60, 255))
    sym_r = QgsSymbol.defaultSymbol(wl.geometryType())
    fl_r  = QgsSimpleFillSymbolLayer()
    fl_r.setColor(ref_fill)
    fl_r.setStrokeColor(border.darker(150))
    fl_r.setStrokeWidth(2.2)
    sym_r.changeSymbolLayer(0, fl_r)

    cats = [
        QgsRendererCategory(0, sym_n, "Secteur réseau"),
        QgsRendererCategory(1, sym_r, "Secteur REF ★"),
    ]
    wl.setRenderer(QgsCategorizedSymbolRenderer("IS_REF", cats))

    QgsProject.instance().addMapLayer(wl)
    n_ref = sum(1 for cid in cell_map if cid == ref_id)
    log.info(f"  → Wedge '{name}': {len(feats)} secteurs  "
             f"(rayon={radius_m}m)  REF={ref_id}")
    return wl


# ═══════════════════════════════════════════════════════════════════════════
#  EXPORT RAPPORT
# ═══════════════════════════════════════════════════════════════════════════

def export_report(P, ref4g, ref3g, ref2g,
                  rows4g, rows_sho, rows_ifho, rows2g):
    """Génère un rapport texte complet dans le dossier projet."""
    now      = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    filename = P + f"rapport_voisins_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

    def stats(rows, label):
        nm = sum(1 for r in rows if r["status"]=="MISSING")
        nc = sum(1 for r in rows if r["status"]=="CONFIGURED")
        np_= sum(1 for r in rows if r["status"]=="POOR_HSR")
        hsr_vals = [r["hsr"] for r in rows if r["hsr"]>0]
        avg_hsr  = round(sum(hsr_vals)/len(hsr_vals),1) if hsr_vals else 0
        lines = [
            f"  {label}",
            f"    Total voisins   : {len(rows)}",
            f"    Manquants (🔴)  : {nm}",
            f"    Configurés (🟢) : {nc}",
            f"    HSR faible (🟠) : {np_}",
            f"    HSR moyen       : {avg_hsr}%",
        ]
        if rows:
            miss_names = [r["tc"] for r in rows if r["status"]=="MISSING"]
            if miss_names:
                lines.append(f"    Voisins manquants :")
                for n in miss_names:
                    lines.append(f"      - {n}")
        return lines

    lines = [
        "="*65,
        "  RAPPORT ANALYSE VOISINS MANQUANTS",
        f"  Généré le : {now}",
        f"  Seuil HSR : {HSR_THRESH}%",
        "="*65,
        "",
        f"  Cellule REF 4G : {ref4g}",
        *stats(rows4g, "4G HO"),
        "",
        f"  Cellule REF 3G : {ref3g}",
        *stats(rows_sho,  "3G SHO"),
        *stats(rows_ifho, "3G IFHO"),
        "",
        f"  Cellule REF 2G : {ref2g}",
        *stats(rows2g, "2G HO"),
        "",
        "="*65,
        "  LÉGENDE",
        "  🔴 ROUGE  = IS_NEIGH=0 (voisin absent de la config réseau)",
        "  🟢 VERT   = IS_NEIGH=1 ET HSR ≥ 97% (voisin bien configuré)",
        "  🟠 ORANGE = IS_NEIGH=1 MAIS HSR < 97% (voisin à optimiser)",
        "  Épaisseur ligne ∝ trafic HO_ATT",
        "="*65,
    ]

    with open(filename, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    log.info(f"  📄 Rapport exporté : {filename}")
    return filename


# ═══════════════════════════════════════════════════════════════════════════
#  CHARGEMENT DES DONNÉES
# ═══════════════════════════════════════════════════════════════════════════
log.info("\n── Chargement des fichiers ──────────────────────────────")

sho  = load_csv(P+CONFIG["FILE_3G_SHO"],
                "3G_SHO",  ["S_CELL","T_CELL","IS_NEIGH","HO_ATT","SHO_HSR","DIST"])
ifho = load_csv(P+CONFIG["FILE_3G_IFHO"],
                "3G_IFHO", ["S_CELL","T_CELL","IS_NEIGH","HO_ATT","IFHO_HSR","DIST"])
ho4g = load_csv(P+CONFIG["FILE_4G_HO"],
                "4G_HO",   ["S_CELL","T_CELL","IS_NEIGH","HO_ATT","DIST"])
ho2g = load_csv(P+CONFIG["FILE_2G_HO"],
                "2G_HO",   ["S_CELL","T_CELL","IS_NEIGH",
                             "OUT_ATT","OUT_FAIL","IN_ATT","IN_FAIL","DIST"])
ref4g= load_csv(P+CONFIG["FILE_4G_CELL_REF"],
                "4G_CELL_REF", ["CELLNAME","X","Y","AZIMUTH"])
ref3g= load_csv(P+CONFIG["FILE_3G_CELL_REF"],
                "3G_CELL_REF", ["CELLNAME","X","Y","AZIMUTH"])
ref2g= load_csv(P+CONFIG["FILE_2G_CELL_REF"],
                "2G_CELL_REF", ["CELLNAME","X","Y","AZIMUTH"])

miss_f = [n for n,d in [("3G_SHO",sho),("3G_IFHO",ifho),("4G_HO",ho4g),
           ("2G_HO",ho2g),("4G_CELL_REF",ref4g),("3G_CELL_REF",ref3g),
           ("2G_CELL_REF",ref2g)] if d is None]
if miss_f:
    raise SystemExit(f"❌ Fichiers manquants: {miss_f}")

db4g = load_xlsx(find_xlsx(P,CONFIG["FILE_DB4G"]), "DB4G")
db3g = load_xlsx(find_xlsx(P,CONFIG["FILE_DB3G"]), "DB3G")
db2g = load_xlsx(find_xlsx(P,CONFIG["FILE_DB2G"]), "DB2G")

# ── Cartes de cellules ────────────────────────────────────────────────────
log.info("\n── Construction des cartes de cellules ──────────────────")

m4g_ref = make_cell_map(ref4g,"CELLNAME","X","Y","AZIMUTH",{"site":"SITENAME"})
m3g_ref = make_cell_map(ref3g,"CELLNAME","X","Y","AZIMUTH",{"site":"SITENAME"})
m2g_ref = make_cell_map(ref2g,"CELLNAME","X","Y","AZIMUTH",{"site":"SITENAME"})

freq_col = "DL FREQ" if "DL FREQ" in ref3g.columns else (
           "DL_FREQ" if "DL_FREQ" in ref3g.columns else None)
if freq_col:
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
            ox,oy=offset_cell(lon,lat,az)
            m4g_new[cid]={"lon":lon,"lat":lat,"az":az,"ox":ox,"oy":oy,
                           "site":str(r["Site_Name"]),"src":"DB"}
        except: pass

m3g_new={}
if db3g is not None:
    for _,r in db3g.iterrows():
        try:
            cid=str(r["UtranCellId"]).strip()
            lon,lat,az=float(r["LONG"]),float(r["LAT"]),float(r["AZIMUTH"])
            ox,oy=offset_cell(lon,lat,az)
            m3g_new[cid]={"lon":lon,"lat":lat,"az":az,"ox":ox,"oy":oy,
                           "site":str(r["3GSiteName"]),
                           "freq":int(r["uarfcnDl"]),"src":"DB"}
        except: pass

m2g_new={}
if db2g is not None:
    sc={}
    for _,r in db2g.iterrows():
        try:
            site=str(r["SiteName"]).strip()
            sec=str(r.get("sector","")).strip()
            cid=f"{site}_{sec}" if sec else f"{site}_S{sc.get(site,0)+1}"
            sc[site]=sc.get(site,0)+1
            lon,lat,az=float(r["LONG"]),float(r["LAT"]),float(r["Azimuth"])
            ox,oy=offset_cell(lon,lat,az)
            m2g_new[cid]={"lon":lon,"lat":lat,"az":az,"ox":ox,"oy":oy,
                           "site":site,"src":"DB"}
        except: pass

m4g_all={**m4g_ref,**m4g_new}
m3g_all={**m3g_ref,**m3g_new}
m2g_all={**m2g_ref,**m2g_new}

log.info(f"  4G: {len(m4g_all)} cellules  |  3G: {len(m3g_all)}  |  2G: {len(m2g_all)}")

# ─────────────────────────────────────────────────────────────────────────────
#  CONSTRUCTION DES FREQ_MAPS  {cell_id: "L800" / "U900" / "G900" / ...}
#  Source prioritaire : CELL_REF  →  complété par DB
# ─────────────────────────────────────────────────────────────────────────────
import re as _re

# ── Correspondances bandes ────────────────────────────────────────────────────
_BAND4G = {'1':'L2100','2':'L1900','3':'L1800','4':'L1700','5':'L850',
           '7':'L2600','8':'L900','12':'L700','17':'L700','20':'L800',
           '28':'L700','38':'TD2600','40':'TD2300','41':'TD2500',
           'B1':'L2100','B3':'L1800','B7':'L2600','B8':'L900',
           'B20':'L800','B28':'L700','B38':'TD2600','B40':'TD2300'}
_BAND3G = {'1':'U2100','5':'U850','8':'U900',
           'B1':'U2100','B5':'U850','B8':'U900'}
_BAND2G = {'GSM900':'G900','GSM1800':'G1800','GSM850':'G850',
           'DCS1800':'G1800','PCS1900':'G1900'}
_UARFCN_U900_DL  = range(2937, 3089)   # Band 8
_UARFCN_U2100_DL = range(10562, 10839) # Band 1
_UARFCN_U850_DL  = range(4357,  4459)  # Band 5

def _uarfcn_to_band(u):
    u = int(u)
    if u in _UARFCN_U900_DL:  return 'U900'
    if u in _UARFCN_U2100_DL: return 'U2100'
    if u in _UARFCN_U850_DL:  return 'U850'
    return f'U({u})'

# ── 4G freq map ───────────────────────────────────────────────────────────────
freq4g = {}
# 1) CELL_REF → FREQ_BAND  (B3→L1800, B7→L2600, …)
for _, _r in ref4g.iterrows():
    _c  = str(_r["CELLNAME"]).strip()
    _fb = str(_r.get("FREQ_BAND","")).strip()
    freq4g[_c] = _BAND4G.get(_fb, _fb if _fb else "?")
# 2) DB4G → freqBand number (integer) ou regex sur nom
if db4g is not None:
    for _, _r in db4g.iterrows():
        _c = str(_r["EUtranCellFDDId"]).strip()
        if _c not in freq4g:
            _fb = str(_r.get("freqBand","")).strip()
            _m  = _re.search(r'(L\d+)', _c, _re.I)
            if _m:
                freq4g[_c] = _m.group(1).upper()
            else:
                freq4g[_c] = _BAND4G.get(_fb, f"B{_fb}" if _fb else "?")

# ── 3G freq map ───────────────────────────────────────────────────────────────
# ⚠️  FREQ_BAND = "B1" pour TOUTES les cellules dans ce fichier (données incorrectes)
#     On utilise DL FREQ (UARFCN) pour distinguer U900 vs U2100
freq3g = {}
_dlfreq_col = "DL FREQ" if "DL FREQ" in ref3g.columns else "DL_FREQ"
for _, _r in ref3g.iterrows():
    _c = str(_r["CELLNAME"]).strip()
    try:
        freq3g[_c] = _uarfcn_to_band(int(_r[_dlfreq_col]))
    except Exception:
        _fb = str(_r.get("FREQ_BAND","")).strip()
        freq3g[_c] = _BAND3G.get(_fb, _fb if _fb else "?")
# DB3G → regex sur nom (contient U900 / U2100) ou uarfcnDl
if db3g is not None:
    for _, _r in db3g.iterrows():
        _c = str(_r["UtranCellId"]).strip()
        if _c not in freq3g:
            _m = _re.search(r'(U\d+)', _c, _re.I)
            if _m:
                freq3g[_c] = _m.group(1).upper()
            else:
                try: freq3g[_c] = _uarfcn_to_band(_r["uarfcnDl"])
                except: freq3g[_c] = "?"

# ── 2G freq map ───────────────────────────────────────────────────────────────
freq2g = {}
for _, _r in ref2g.iterrows():
    _c  = str(_r["CELLNAME"]).strip()
    _fb = str(_r.get("FREQ_BAND","")).strip()
    freq2g[_c] = _BAND2G.get(_fb, _fb if _fb else "?")
if db2g is not None:
    _sc = {}
    for _, _r in db2g.iterrows():
        _site = str(_r["SiteName"]).strip()
        _sec  = str(_r.get("sector","")).strip()
        _cid  = f"{_site}_{_sec}" if _sec else f"{_site}_S{_sc.get(_site,0)+1}"
        _sc[_site] = _sc.get(_site, 0) + 1
        if _cid not in freq2g:
            _sys = str(_r.get("cSysType","")).strip()
            freq2g[_cid] = _BAND2G.get(_sys, _sys if _sys else "G?")

log.info(f"  Freq4G: {len(freq4g)} entrées  |  "
         f"Freq3G: {len(freq3g)}  |  Freq2G: {len(freq2g)}")
log.info(f"  4G bandes: {sorted(set(freq4g.values()))}")
log.info(f"  3G bandes: {sorted(set(freq3g.values()))}")
log.info(f"  2G bandes: {sorted(set(freq2g.values()))}")

# ─────────────────────────────────────────────────────────────────────────────
#  Listes pour menus  — format : "NOM_CELLULE  [L800]"
# ─────────────────────────────────────────────────────────────────────────────
def _tag(cell, fmap):
    """Retourne 'CELLNAME  [BAND]' ou 'CELLNAME  [?]' si inconnu."""
    b = fmap.get(cell.strip(), "?")
    return f"{cell}  [{b}]"

_raw_4g = sorted(ho4g["S_CELL"].dropna().unique().tolist())
_raw_3g = sorted(sho ["S_CELL"].dropna().unique().tolist())
_raw_2g = sorted(ho2g["S_CELL"].dropna().unique().tolist())

cells_4g         = [_tag(c, freq4g) for c in _raw_4g]
cells_3g         = [_tag(c, freq3g) for c in _raw_3g]
cells_2g         = [_tag(c, freq2g) for c in _raw_2g]
# Maps display → original  (pour récupérer le vrai nom après sélection)
_disp2raw_4g = {_tag(c,freq4g): c for c in _raw_4g}
_disp2raw_3g = {_tag(c,freq3g): c for c in _raw_3g}
_disp2raw_2g = {_tag(c,freq2g): c for c in _raw_2g}


# ═══════════════════════════════════════════════════════════════════════════
#  FENÊTRE DE SÉLECTION
# ═══════════════════════════════════════════════════════════════════════════
class CellSelector(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Analyse voisins manquants — Sélection des cellules REF")
        self.setMinimumWidth(520)
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)

        lay = QVBoxLayout(self)
        lay.setSpacing(10)
        lay.setContentsMargins(16,16,16,16)

        # Titre
        title = QLabel("📡  Choisissez une cellule de référence par technologie  —  [bande fréquence]")
        title.setStyleSheet(
            "font-size:14px;font-weight:bold;padding:8px;"
            "background:#f0f4ff;border-radius:4px;")
        title.setAlignment(Qt.AlignCenter)
        lay.addWidget(title)

        # Info seuil
        info = QLabel(
            f"  Seuil HSR configuré : {HSR_THRESH}%\n"
            f"  🟢 Vert ≥ {HSR_THRESH}%   |   "
            f"🟠 Orange < {HSR_THRESH}%   |   "
            f"🔴 Rouge = manquant (IS_NEIGH=0)")
        info.setStyleSheet(
            "color:#333;font-size:11px;padding:6px 10px;"
            "background:#fffbe6;border:1px solid #e0c060;border-radius:4px;")
        lay.addWidget(info)

        def make_group(title_txt, color, cells, attr_s, attr_c):
            grp = QGroupBox(f"  {title_txt}")
            grp.setStyleSheet(
                f"QGroupBox{{font-weight:bold;color:{color};"
                f"font-size:12px;border:1.5px solid {color};"
                f"border-radius:5px;margin-top:8px;padding:4px;}}"
                f"QGroupBox::title{{subcontrol-origin:margin;"
                f"left:10px;padding:0 4px;}}")
            v = QVBoxLayout(grp)
            search = QLineEdit()
            search.setPlaceholderText("🔍  Filtrer…")
            search.setStyleSheet("padding:4px;font-size:11px;")
            combo  = QComboBox()
            combo.setMaxVisibleItems(18)
            combo.setStyleSheet("font-size:11px;padding:2px;")
            combo.addItems(cells)
            search.textChanged.connect(
                lambda t,c=combo,cl=cells: (
                    c.clear(),
                    c.addItems([x for x in cl if t.lower() in x.lower()])
                ))
            v.addWidget(search); v.addWidget(combo)
            # Tooltip explaining the [BAND] tag
            combo.setToolTip(
                "Format affiché : NOM_CELLULE  [BANDE]\n"
                "Filtrez par nom ou par bande (ex: L800, U900, G900)"
            )
            setattr(self, attr_s, search)
            setattr(self, attr_c, combo)
            return grp

        lay.addWidget(make_group("4G — Cellule de référence  (ex: L800 · L1800 · L2100 · L2600)",
                                 "#0055cc", cells_4g, "s4","c4"))
        lay.addWidget(make_group("3G — Cellule de référence  (ex: U900 · U2100)  —  SHO + IFHO",
                                 "#cc5500", cells_3g, "s3","c3"))
        lay.addWidget(make_group("2G — Cellule de référence  (ex: G900 · G1800)",
                                 "#006633", cells_2g, "s2","c2"))

        # Boutons
        btn_row = QHBoxLayout()
        btn_ok  = QPushButton("▶   Lancer l'analyse")
        btn_ok.setMinimumHeight(38)
        btn_ok.setStyleSheet(
            "QPushButton{background:#0055cc;color:white;font-weight:bold;"
            "padding:8px 24px;border-radius:5px;font-size:13px;}"
            "QPushButton:hover{background:#003fa3;}"
            "QPushButton:pressed{background:#002d80;}")
        btn_ok.clicked.connect(self.accept)
        btn_cancel = QPushButton("Annuler")
        btn_cancel.setMinimumHeight(38)
        btn_cancel.setStyleSheet("padding:8px 16px;font-size:12px;"
                                 "border-radius:5px;")
        btn_cancel.clicked.connect(self.reject)
        btn_row.addStretch()
        btn_row.addWidget(btn_ok)
        btn_row.addWidget(btn_cancel)
        lay.addLayout(btn_row)

    def cells(self):
        """Retourne les vrais noms de cellules (sans le tag [BAND])."""
        disp4 = self.c4.currentText()
        disp3 = self.c3.currentText()
        disp2 = self.c2.currentText()
        raw4  = _disp2raw_4g.get(disp4, disp4.split("  [")[0])
        raw3  = _disp2raw_3g.get(disp3, disp3.split("  [")[0])
        raw2  = _disp2raw_2g.get(disp2, disp2.split("  [")[0])
        return raw4, raw3, raw2


dlg = CellSelector()
if not dlg.exec_():
    raise SystemExit("Annulé.")

REF_4G, REF_3G, REF_2G = dlg.cells()
log.info(f"\n  4G → {REF_4G}  [{freq4g.get(REF_4G,'?')}]")
log.info(f"  3G → {REF_3G}  [{freq3g.get(REF_3G,'?')}]")
log.info(f"  2G → {REF_2G}  [{freq2g.get(REF_2G,'?')}]\n")


# ═══════════════════════════════════════════════════════════════════════════
#  BARRE DE PROGRESSION
# ═══════════════════════════════════════════════════════════════════════════
TOTAL_STEPS = 18
progress = QProgressDialog("Analyse en cours…", "Annuler", 0, TOTAL_STEPS)
progress.setWindowTitle("Analyse voisins manquants")
progress.setWindowModality(Qt.WindowModal)
progress.setMinimumWidth(420)
progress.setMinimumDuration(0)
progress.setAutoClose(False)
progress.setValue(0)

def step(msg, n=1):
    progress.setValue(progress.value() + n)
    progress.setLabelText(msg)
    QApplication.processEvents()
    if progress.wasCanceled():
        raise SystemExit("Annulé par l'utilisateur.")


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSE 4G
# ═══════════════════════════════════════════════════════════════════════════
step(f"⚙️  4G — calcul HO HSR…")
log.info("─"*65)
log.info(f"  4G : {REF_4G}")
g4_rows = []

if REF_4G in m4g_all:
    si = m4g_all[REF_4G]
    if "HO_FAIL" in ho4g.columns:
        ho4g["HO_HSR"] = ho4g.apply(lambda r: hsr_safe(r["HO_ATT"],r["HO_FAIL"]),axis=1)
    elif "HO_HSR" not in ho4g.columns:
        ho4g["HO_HSR"] = 0.0

    known=set()
    for _,r in ho4g[ho4g["S_CELL"]==REF_4G].iterrows():
        tc=r["T_CELL"]; isn=int(r["IS_NEIGH"])
        att=int(r["HO_ATT"]); hsr=float(r.get("HO_HSR",0)); dist=float(r["DIST"])
        ti=m4g_all.get(tc)
        if ti: tox,toy=ti["ox"],ti["oy"]
        else:
            try: tox,toy=offset_cell(float(r["T_X"]),float(r["T_Y"]),0)
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
    log.info(f"  4G total: {len(g4_rows)}  "
             f"[🔴{sum(1 for r in g4_rows if r['status']=='MISSING')} "
             f"| 🟢{sum(1 for r in g4_rows if r['status']=='CONFIGURED')} "
             f"| 🟠{sum(1 for r in g4_rows if r['status']=='POOR_HSR')}]")
else:
    log.error(f"  ✗ '{REF_4G}' introuvable dans les données 4G")

step("🗺️  4G — création des couches…")
add_wedge_layer("Wedges_4G", REF_4G, m4g_all, QColor(30,160,230,100), RADIUS_4G)
step("🗺️  4G — points cellules…")
add_all_cell_points("Cell_Points_4G", m4g_all, QColor(0,120,220), "4G", REF_4G,
                    freq_map=freq4g)
step("🗺️  4G — lignes voisins…")
add_neighbor_lines("REF_Neighbors_4G", g4_rows)


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSE 3G SHO
# ═══════════════════════════════════════════════════════════════════════════
step(f"⚙️  3G SHO — analyse…")
log.info("─"*65)
log.info(f"  3G : {REF_3G}")
sho_rows=[]; ref_freq=0; si3=None

if REF_3G in m3g_all:
    si3=m3g_all[REF_3G]; ref_freq=si3.get("freq",0)
    known=set()
    for _,r in sho[sho["S_CELL"]==REF_3G].iterrows():
        tc=r["T_CELL"]; isn=int(r["IS_NEIGH"])
        att=int(r["HO_ATT"]); hsr=float(r.get("SHO_HSR",0)); dist=float(r["DIST"])
        ti=m3g_all.get(tc)
        if ti: tox,toy=ti["ox"],ti["oy"]
        else:
            try: tox,toy=offset_cell(float(r["T_X"]),float(r["T_Y"]),0)
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
    log.info(f"  SHO: {len(sho_rows)}  "
             f"[🔴{sum(1 for r in sho_rows if r['status']=='MISSING')} "
             f"| 🟢{sum(1 for r in sho_rows if r['status']=='CONFIGURED')}]")
else:
    log.error(f"  ✗ '{REF_3G}' introuvable 3G")

# ── IFHO ─────────────────────────────────────────────────────────────────
step(f"⚙️  3G IFHO — analyse…")
ifho_rows=[]
if si3 is not None:
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
            try: tox,toy=offset_cell(float(r["T_X"]),float(r["T_Y"]),0)
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
    log.info(f"  IFHO: {len(ifho_rows)}  "
             f"[🔴{sum(1 for r in ifho_rows if r['status']=='MISSING')} "
             f"| 🟢{sum(1 for r in ifho_rows if r['status']=='CONFIGURED')}]")

step("🗺️  3G — création des couches…")
add_wedge_layer("Wedges_3G", REF_3G, m3g_all, QColor(140,140,140,100), RADIUS_3G)
step("🗺️  3G — points cellules…")
add_all_cell_points("Cell_Points_3G", m3g_all, QColor(220,100,0), "3G", REF_3G,
                    freq_map=freq3g)
step("🗺️  3G — lignes SHO…")
add_neighbor_lines("REF_Neighbors_3G_SHO",  sho_rows)
step("🗺️  3G — lignes IFHO…")
add_neighbor_lines("REF_Neighbors_3G_IFHO", ifho_rows)


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSE 2G
# ═══════════════════════════════════════════════════════════════════════════
step(f"⚙️  2G — analyse…")
log.info("─"*65)
log.info(f"  2G : {REF_2G}")
g2_rows=[]; si2=None

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
            try: tox,toy=offset_cell(float(r["T_X"]),float(r["T_Y"]),0)
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
    log.info(f"  2G total: {len(g2_rows)}  "
             f"[🔴{sum(1 for r in g2_rows if r['status']=='MISSING')} "
             f"| 🟢{sum(1 for r in g2_rows if r['status']=='CONFIGURED')} "
             f"| 🟠{sum(1 for r in g2_rows if r['status']=='POOR_HSR')}]")
else:
    log.error(f"  ✗ '{REF_2G}' introuvable 2G")

step("🗺️  2G — création des couches…")
add_wedge_layer("Wedges_2G", REF_2G, m2g_all, QColor(220,80,160,100), RADIUS_2G)
step("🗺️  2G — points cellules…")
add_all_cell_points("Cell_Points_2G", m2g_all, QColor(0,170,60), "2G", REF_2G,
                    freq_map=freq2g)
step("🗺️  2G — lignes voisins…")
add_neighbor_lines("REF_Neighbors_2G", g2_rows)


# ═══════════════════════════════════════════════════════════════════════════
#  ZOOM AUTOMATIQUE + EXPORT RAPPORT + POPUP RÉSUMÉ
# ═══════════════════════════════════════════════════════════════════════════
step("🔍  Zoom sur la cellule de référence…")

# Zoom sur cellule 4G REF (ou 3G si absente)
zoom_done = False
for ref_id, cmap in [(REF_4G, m4g_all),(REF_3G, m3g_all),(REF_2G, m2g_all)]:
    if ref_id in cmap:
        zoom_to_cell(cmap, ref_id, margin_deg=0.025)
        zoom_done = True
        break
if not zoom_done:
    from qgis.utils import iface
    iface.mapCanvas().refresh()

step("📄  Export du rapport…")
rapport_file = export_report(P, REF_4G, REF_3G, REF_2G,
                              g4_rows, sho_rows, ifho_rows, g2_rows)

step("✅  Terminé !")
progress.close()

log.info("\n" + "="*65)
log.info("  ✅  ANALYSE TERMINÉE")
log.info("="*65)

# ── Popup résumé final ────────────────────────────────────────────────────
def s(rows, key):
    return sum(1 for r in rows if r["status"]==key)

hsr_vals_4g  = [r["hsr"] for r in g4_rows   if r["hsr"]>0]
hsr_vals_sho = [r["hsr"] for r in sho_rows  if r["hsr"]>0]
hsr_vals_2g  = [r["hsr"] for r in g2_rows   if r["hsr"]>0]
avg4  = round(sum(hsr_vals_4g )/len(hsr_vals_4g ),1) if hsr_vals_4g  else 0
avg3  = round(sum(hsr_vals_sho)/len(hsr_vals_sho),1) if hsr_vals_sho else 0
avg2  = round(sum(hsr_vals_2g )/len(hsr_vals_2g ),1) if hsr_vals_2g  else 0

summary = f"""
╔══════════════════════════════════════════════════╗
║     RÉSULTATS ANALYSE VOISINS MANQUANTS          ║
╠══════════════════════════════════════════════════╣
║                                                  ║
║  4G ▶ {REF_4G[:38]:<38}  ║
║       Bande : {freq4g.get(REF_4G,'?'):<37}  ║
║       Voisins : {len(g4_rows):<4}  🔴 {s(g4_rows,'MISSING'):<3} manquants
║       🟢 {s(g4_rows,'CONFIGURED'):<3} configurés  🟠 {s(g4_rows,'POOR_HSR'):<3} HSR faible
║       HSR moyen : {avg4}%
║                                                  ║
║  3G ▶ {REF_3G[:38]:<38}  ║
║       Bande : {freq3g.get(REF_3G,'?'):<37}  ║
║       SHO  : {len(sho_rows):<4}  🔴 {s(sho_rows,'MISSING'):<3} manquants
║       IFHO : {len(ifho_rows):<4}  🔴 {s(ifho_rows,'MISSING'):<3} manquants
║       HSR SHO moyen : {avg3}%
║                                                  ║
║  2G ▶ {REF_2G[:38]:<38}  ║
║       Bande : {freq2g.get(REF_2G,'?'):<37}  ║
║       Voisins : {len(g2_rows):<4}  🔴 {s(g2_rows,'MISSING'):<3} manquants
║       🟢 {s(g2_rows,'CONFIGURED'):<3} configurés  🟠 {s(g2_rows,'POOR_HSR'):<3} HSR faible
║       HSR moyen : {avg2}%
║                                                  ║
║  📄 Rapport : {os.path.basename(rapport_file):<34}  ║
║  📋 Log     : {os.path.basename(LOG_FILE):<34}  ║
╚══════════════════════════════════════════════════╝

LÉGENDE DES COUCHES :
  🔴 ROUGE  = voisin MANQUANT  (IS_NEIGH=0)
  🟢 VERT   = voisin CONFIGURÉ (HSR ≥ {HSR_THRESH}%)
  🟠 ORANGE = voisin en place  (HSR < {HSR_THRESH}%)
  ★  Étoile = cellule de référence
  Épaisseur ligne ∝ trafic HO_ATT
"""

msg = QMessageBox()
msg.setWindowTitle("✅  Analyse terminée — PFE Tunisie Telecom 2026")
msg.setText(summary)
msg.setIcon(QMessageBox.Information)
msg.setStandardButtons(QMessageBox.Ok)
msg.setMinimumWidth(560)
msg.exec_()

log.info(f"  Rapport : {rapport_file}")
log.info(f"  Log     : {LOG_FILE}")
