import ngsolve as ngs
import netgen.occ as occ
from netgen.webgui import Draw as DrawGeo


def Mesh(occ_geometry,
         maxh = 1.,
         *args, **kwargs) -> ngs.Mesh:
    "Generate a NGSolve mesh from a Netgen OCC geometry"
    mesh = ngs.Mesh(occ_geometry.GenerateMesh(maxh=maxh, *args, **kwargs))
    return mesh


def square(side : float = 1.,
           inclusion : bool = False,
           draw : bool = False) -> occ.OCCGeometry:
    """
    Create an OpenCASCADE (OCC) geometry of a square.
    """
    box = occ.Rectangle(side, side).Face()
    box.name = "out"
    if inclusion:
        box_in = occ.MoveTo(side/4,side/4).Rectangle(side/2, side/2).Face()
        box_in.name = "in"
        geo = occ.Glue([box-box_in, box_in])
    else:
        geo = box
    if draw:
        DrawGeo(geo)
    return occ.OCCGeometry(geo, dim=2)

def disk(radius: float = 1.,
         inclusion : bool = False,
         draw: bool = False) -> occ.OCCGeometry:
    """
    Create an OpenCASCADE (OCC) geometry of a disk.
    """
    box = occ.Circle(occ.Pnt(0, 0), radius).Face()
    box.name = "out"
    if inclusion:
        box_in = occ.Circle(occ.Pnt(0, 0), radius/2).Face()
        box_in.name = "in"
        geo = occ.Glue([box-box_in, box_in])
    else:
        geo = box
    if draw:
        DrawGeo(geo)
    return occ.OCCGeometry(geo, dim=2)


### Provisoire

from ngsolve import Mesh, unit_square
from netgen.geom2d import CSG2d, Circle, Rectangle
from netgen.geom2d import EdgeInfo as EI, PointInfo as PI, Solid2d


def gaped_inductor(
    airgap: float = 1e-3,
    h: float = 0.005,
) -> Mesh:
    """
    Create a mesh of an inductor core with (optionally) an air gap.

    Uses geom2d's Rectangle with pmin/pmax coordinates, boolean subtraction,
    and CSG2d geometry container.

    The geometry consists of:
      - A bounding "air" region
      - An "iron" core region
      - An air cavity inside the iron ("airIron")
      - Two conductors (positive and negative)
      - Optionally, a thin air gap between iron sections

    Parameters
    ----------
    airgap : float, optional
        Width of the air gap (default 1 mm). Set to 0 to remove it.
    h : float, optional
        Maximum element size in the air region (default 5 mm).

    Returns
    -------
    Mesh
        An NGSolve mesh of the inductor geometry.
    """
    # ------------------------------------------------------------------
    # 1. Define 2D primitives using geom2d rectangles
    #    Rectangle(pmin=(x1,y1), pmax=(x2,y2), mat=..., bc=...)
    #    - mat  : material name
    #    - bc   : boundary condition name (applied to specified edges)
    # ------------------------------------------------------------------
    box_ = Rectangle(
        pmin=(0.02, 0.02), pmax=(0.08, 0.08),
        mat="air", bc="out",
    )

    iron_ = Rectangle(
        pmin=(0.04, 0.04), pmax=(0.06, 0.06),
        mat="iron",
    )

    airIron_ = Rectangle(
        pmin=(0.045, 0.045), pmax=(0.055, 0.055),
        mat="air",
    )

    positiveCond_ = Rectangle(
        pmin=(0.05, 0.045), pmax=(0.055, 0.055),
        mat="condP",
    )

    negativeCond_ = Rectangle(
        pmin=(0.06, 0.045), pmax=(0.065, 0.055),
        mat="condN",
    )

    # ------------------------------------------------------------------
    # 2. Build composite geometry using boolean subtraction (-)
    #    geom2d Solid2d supports intuitive +, -, & operators for CSG.
    # ------------------------------------------------------------------
    if airgap > 0:
        # Create the air gap rectangle: thin slice at x = [0.03, 0.046], y centred at 0.05
        airgap_ = Rectangle(
            pmin=(0.03, 0.05 - airgap / 2),
            pmax=(0.046, 0.05 + airgap / 2),
            mat="air",
        )
        # Iron core = full iron block minus its air cavity, gap, and conductors
        iron = iron_ - airIron_ - airgap_ - negativeCond_ - positiveCond_
        # Air domain = outer box minus iron and conductors
        air = box_ - iron - positiveCond_ - negativeCond_
    else:
        # No air gap: simpler subtraction
        iron = iron_ - airIron_ - negativeCond_ - positiveCond_
        air = box_ - iron - negativeCond_ - positiveCond_

    # ------------------------------------------------------------------
    # 3. Set material properties and mesh sizing
    #    .Mat(name)         — set material name (redundant here but explicit)
    #    .Maxh(value)       — set maximum element size for this region
    # ------------------------------------------------------------------
    air.Mat("air")
    iron.Mat("iron").Maxh(h / 10)  # Finer mesh in iron

    # ------------------------------------------------------------------
    # 4. Add all domains to a CSG2d geometry container and generate the mesh
    # ------------------------------------------------------------------
    geo = CSG2d()
    geo.Add(air)
    geo.Add(iron)
    geo.Add(positiveCond_)
    geo.Add(negativeCond_)

    return Mesh(geo.GenerateMesh(maxh=h))
