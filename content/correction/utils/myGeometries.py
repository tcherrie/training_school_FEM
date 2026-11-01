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
