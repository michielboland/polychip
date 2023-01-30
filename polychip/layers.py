import functools
import shapely
import shapely.geometry
import shapely.wkt

from enum import Enum, unique
from svg_parse import *

def coerce_multipoly(poly):
    """Coerces the given poly to be a MultiPoly.

    Results from shapely operations can be Polygons or MultiPolygons. It's convenient to
    have them always be MultiPolygons so that, for example, we can iterate over whatever it
    is that was returned.

    Args:
        poly (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): The poly to coerce.

    Returns:
        shapely.geometry.MultiPolygon: The coerced poly.
    """
    if type(poly) == shapely.geometry.Polygon:
        return shapely.geometry.MultiPolygon([poly])
    return poly


@unique
class Layer(Enum):
    """ Types of layers."""
    METAL = "Metal"
    POLY = "Poly"
    DIFF = "Diff"
    CONTACTS = "Contacts"
    QNAMES = "QNames"
    SNAMES = "SNames"
    PNAMES = "PNames"
    CAPACITORS = "Capacitors"

    def path(self):
        return "./svg:g[@inkscape:groupmode='layer'][@inkscape:label='" + self.value + "']"


class Label(object):
    """Represents the text and extents of a label.

    Args:
    Attributes:
        text (str): The text of the label.
        extents (shapely.geometry.LineString): A line from bottom beginning
            to top end (relative to the text string, not its orientation).
    """
    def __init__(self, text, extents):
        self.text = text
        self.extents = extents
        self.center = extents.centroid

    def to_dict(self):
        """Converts to a dictionary, for JSON encoding."""
        return {
            "__POLYCHIP_OBJECT__": "Label",
            "text": self.text,
            "extents": self.extents.wkt,
        }

    @staticmethod
    def from_dict(d):
        """Converts a dictionary to a Transistor, for JSON decoding."""
        assert d["__POLYCHIP_OBJECT__"] == "Label", "Label.from_dict wasn't given its expected dict: " + str(d)
        return Label(d["text"], shapely.wkt.loads(d["extents"]))
        

class InkscapeFile:
    """Represents all the paths and names found in an Inkscape file.

    Args:
        root (xml.etree.ElementTree.Element): The root element for the Inkscape document.

    Attributes:
        contacts (shapely.geometry.MultiPoint): All the found contacts. Each point
            represents the center position of a rectangular contact.
        qnames([Label]): The list of found transistor labels.
        snames([Label]): The list of found signal labels.
        pnames([Label]): The list of found pin labels.
        poly_array([shapely.geometry.Polygon]): The list of found polysilicon polygons.
        metal_array([shapely.geometry.Polygon]): The list of found metal polygons.
        diff_array([shapely.geometry.Polygon]): The list of found diff polygons.
            Note that this will be altered in a later stage when transistors are found.
        multicontact(shapely.geometry.MultiPolygon): All the found contact polygons.
        multipoly(shapely.geometry.MultiPolygon): All the found polysilicon polygons.
        multidiff(shapely.geometry.MultiPolygon): All the found diffusion polygons.
        multimetal(shapely.geometry.MultiPolygon): All the found metal polygons.
        multicaps(shapely.geometry.MultiPolygon): All the found capacitor polygons.
    """
    def __init__(self, root):
        self.qnames = []
        self.snames = []
        self.pnames = []
        self.contact_array = []
        self.poly_array = []
        self.metal_array = []
        self.diff_array = []
        self.multicontact = shapely.geometry.MultiPolygon()
        self.multipoly = shapely.geometry.MultiPolygon()
        self.multidiff = shapely.geometry.MultiPolygon()
        self.multimetal = shapely.geometry.MultiPolygon()
        self.multicaps = shapely.geometry.MultiPolygon()

        self.to_screen_coords_transform_ = self.extract_screen_transform(root)

        self.transform = {l: self.to_screen_coords_transform_ for l in Layer}

        self.contact_paths = {}
        poly_paths = {}
        diff_paths = {}
        metal_paths = {}
        capacitor_paths = {}

        layer = {}

        for l in (Layer.POLY, Layer.DIFF, Layer.METAL, Layer.CONTACTS):
            layer[l] = root.findall(l.path(), namespaces)[0]

        for l in (Layer.QNAMES, Layer.SNAMES, Layer.PNAMES, Layer.CAPACITORS):
            namelayer = root.findall(l.path(), namespaces)
            if len(namelayer) > 0:
                layer[l] = namelayer[0]

        for y in (y for y in Layer if y in layer):
            t = Transform.parse(layer[y].get('transform'))
            self.transform[y] = self.transform[y] @ t
        shapes = {}

        for l in (
            Layer.POLY,
            Layer.DIFF,
            Layer.METAL,
            Layer.CONTACTS,
            Layer.CAPACITORS,
        ):
            shapes[l] = root.findall(l.path() + "/svg:path", namespaces)
            shapes[l] += root.findall(l.path() + "/svg:rect", namespaces)

        for l in (Layer.QNAMES, Layer.SNAMES, Layer.PNAMES):
            shapes[l] = root.findall(l.path() + "/svg:text", namespaces)

        print("Processing {:d} contact paths".format(len(shapes[Layer.CONTACTS])))
        for p in shapes[Layer.CONTACTS]:
            self.contact_paths['c_' + p.get('id')] = svgelement_to_shapely_polygon(p, self.transform[Layer.CONTACTS])

        print("Processing {:d} poly paths".format(len(shapes[Layer.POLY])))
        for p in shapes[Layer.POLY]:
            poly_paths['p_' + p.get('id')] = svgelement_to_shapely_polygon(p, self.transform[Layer.POLY])

        print("Processing {:d} diff paths".format(len(shapes[Layer.DIFF])))
        for p in shapes[Layer.DIFF]:
            diff_paths['p_' + p.get('id')] = svgelement_to_shapely_polygon(p, self.transform[Layer.DIFF])

        print("Processing {:d} metal paths".format(len(shapes[Layer.METAL])))
        for p in shapes[Layer.METAL]:
            metal_paths['p_' + p.get('id')] = svgelement_to_shapely_polygon(p, self.transform[Layer.METAL])

        print("Processing {:d} capacitor paths".format(len(shapes[Layer.CAPACITORS])))
        for p in shapes[Layer.CAPACITORS]:
            capacitor_paths['p_' + p.get('id')] = svgelement_to_shapely_polygon(p, self.transform[Layer.CAPACITORS])

        print("Processing qnames text")
        for t in shapes[Layer.QNAMES]:
            text, extents = parse_shapely_text(t, self.transform[Layer.QNAMES])
            self.qnames.append(Label(text, extents))

        print("Processing snames text")
        for t in shapes[Layer.SNAMES]:
            text, extents = parse_shapely_text(t, self.transform[Layer.SNAMES])
            self.snames.append(Label(text, extents))

        print("Processing pnames text")
        for t in shapes[Layer.PNAMES]:
            text, extents = parse_shapely_text(t, self.transform[Layer.PNAMES])
            self.pnames.append(Label(text, extents))
            self.snames.append(Label(text, extents))

        print("Merging overlapping sections. Before merge:")
        print("{:d} contacts".format(len(self.contact_paths)))
        print("{:d} diffs".format(len(diff_paths)))
        print("{:d} polys".format(len(poly_paths)))
        print("{:d} metals".format(len(metal_paths)))
        print("After merging:")
        # print(diff_paths['p_rect10018'])

        self.multicontact = coerce_multipoly(shapely.ops.unary_union(
            [p for p in self.contact_paths.values() if p is not None]))
        self.contact_array = list(self.multicontact.geoms)
        list.sort(self.contact_array, key = functools.cmp_to_key(InkscapeFile.poly_cmp))
        print("{:d} contacts".format(len(self.contact_array)))

        self.multidiff = coerce_multipoly(shapely.ops.unary_union(
            [p for p in diff_paths.values() if p is not None]))
        self.diff_array = list(self.multidiff.geoms)
        list.sort(self.diff_array, key = functools.cmp_to_key(InkscapeFile.poly_cmp))
        print("{:d} diffs".format(len(self.diff_array)))

        self.multipoly = coerce_multipoly(shapely.ops.unary_union(
            [p for p in poly_paths.values() if p is not None]))
        self.poly_array = list(self.multipoly.geoms)
        list.sort(self.poly_array, key = functools.cmp_to_key(InkscapeFile.poly_cmp))
        print("{:d} polys".format(len(self.poly_array)))

        self.multicaps = coerce_multipoly(shapely.ops.unary_union(
            [p for p in capacitor_paths.values() if p is not None]))

        self.multimetal = coerce_multipoly(shapely.ops.unary_union(
            [p for p in metal_paths.values() if p is not None]))
        self.metal_array = list(self.multimetal.geoms)
        list.sort(self.metal_array, key = functools.cmp_to_key(InkscapeFile.poly_cmp))
        print("{:d} metals".format(len(self.metal_array)))

        print("{:d} qnames".format(len(self.qnames)))
        print("{:d} snames".format(len(self.snames)))
        print("{:d} pnames".format(len(self.pnames)))


    def extract_screen_transform(self, root):
        """Extracts the height, in pixels, of the document.

        Args:
            root (xml.etree.ElementTree.Element): The root element for the Inkscape document.

        Returns:
            Transform: The transform to get from SVG coordinates to Inkscape screen coordinates.
        """
        height = root.get('height')
        width = root.get('width')
        xdpi = root.get(qname(root, "inkscape:export-xdpi"))
        ydpi = root.get(qname(root, "inkscape:export-ydpi"))
        if height.endswith('mm'):
            h = float(xdpi) * float(height[:-2]) / 25.4
        else:
            h = float(height)

        if width.endswith('mm'):
            w = float(ydpi) * float(width[:-2]) / 25.4
        else:
            w = float(width)

        transform = Transform(1, 0, 0, -1, 0, h)

        # If there's a viewBox, then the document scale needs adjustment.
        viewbox = root.get('viewBox')
        if viewbox is None:
            return transform
        extents = [float(x) for x in re.split('[, ]', viewbox)]
        scalex = w / extents[2]
        scaley = h / extents[3]
        return transform @ Transform.scale(scalex, scaley)


    def replace_diff_array(self, diffs):
        """Replaces the drawing's diff_array with the given one, sorting it first.

        This happens after transistors are identified. The transistor gates split existing
        diffs in two.

        Args:
            diffs ([shapely.geometry.Polygon]): The array of diff polygons.
        """
        self.diff_array = diffs
        list.sort(self.diff_array, key = functools.cmp_to_key(InkscapeFile.poly_cmp))


    @staticmethod
    def poly_cmp(poly1, poly2):
        """Provides an ordering for two polygons based on their bounding box.

        The polygon whose bounding box is leftmost of the two is the lower one. If both polygons are
        left-aligned, then the polygon that is lowermost of the two is the lower one.

        Args:
            poly1 (shapely.geometry.Polygon): The first polygon to compare.
            poly2 (shapely.geometry.Polygon): The polygon to compare the first polygon to.

        Returns:
            int:
                -1 if poly1 is "less than" poly2
                1 if poly1 is "greater than" poly2
                0 if poly1 is "equal to" poly 2

        """
        minx1, miny1, _, _ = poly1.bounds
        minx2, miny2, _, _ = poly2.bounds
        if minx1 < minx2:
            return -1
        if minx1 > minx2:
            return 1
        if miny1 < miny2:
            return -1
        if miny1 > miny2:
            return 1
        return 0
