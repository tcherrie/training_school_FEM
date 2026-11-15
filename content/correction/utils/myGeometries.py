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

from netgen.geom2d import CSG2d, Rectangle
from time import time

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
    iron.Mat("iron")

    # ------------------------------------------------------------------
    # 4. Add all domains to a CSG2d geometry container and generate the mesh
    # ------------------------------------------------------------------
    geo = CSG2d()
    geo.Add(air)
    geo.Add(iron)
    geo.Add(positiveCond_)
    geo.Add(negativeCond_)
    
    t0 = time()
    mesh = ngs.Mesh(geo.GenerateMesh(maxh=h))
    print(f"Mesh generated in {time()-t0 :.3f} s ({mesh.ne} elements)")
    return mesh


######################
# 3D 
"""
Created on Mon Oct 23 21:23:29 2023

@author: A.Cesarano, T.Cherriere, T. Gauthey, M. Hage-Hassan
"""


#from netgen.csg import *
#from netgen.occ import *
from netgen.occ import Pnt, Segment, ArcOfCircle, Segment, Face, Wire, X, Y, Z, Glue, OCCGeometry, Axis
from ngsolve import sqrt
#from ngsolve.internal import *
from netgen.webgui import Draw as DrawGeo

outer_len = 150*10**(-3)
coil_z_len = 40*10**(-3)
coil_rad_outer = 30*10**(-3)
coil_rad_inner = 20*10**(-3)
core_d = 5*10**(-3)
core_air_z = 30*10**(-3)
core_air_x = 30*10**(-3)

splitting = 8


def geometryCoil(splitting = 1, outer_len = 150e-3, coil_z_len = 46e-3, 
                 coil_rad_outer = 24e-3, coil_rad_inner = 20e-3, core_d = 5e-3, 
                 core_air_z = 25e-3, core_air_x = 20e-3):
    
    if splitting == 8:
        o1 = Pnt( 0.5*outer_len,  0.5*outer_len, 0)
        o2 = Pnt( 0            ,  0.5*outer_len, 0)
        o3 = Pnt( 0            ,  0            , 0)
        o4 = Pnt( 0.5*outer_len,  0            , 0)
    elif splitting == 4:
        o1 = Pnt( 0.5*outer_len,  0.5*outer_len, -0.5*outer_len)
        o2 = Pnt( 0            ,  0.5*outer_len, -0.5*outer_len)
        o3 = Pnt( 0            ,  0            , -0.5*outer_len)
        o4 = Pnt( 0.5*outer_len,  0            , -0.5*outer_len)
    elif splitting == 2:
        o1 = Pnt( 0.5*outer_len,  0.5*outer_len, -0.5*outer_len)
        o2 = Pnt(-0.5*outer_len,  0.5*outer_len, -0.5*outer_len)
        o3 = Pnt(-0.5*outer_len,  0            , -0.5*outer_len)
        o4 = Pnt( 0.5*outer_len,  0            , -0.5*outer_len)
    else:
        o1 = Pnt( 0.5*outer_len,  0.5*outer_len, -0.5*outer_len)
        o2 = Pnt(-0.5*outer_len,  0.5*outer_len, -0.5*outer_len)
        o3 = Pnt(-0.5*outer_len, -0.5*outer_len, -0.5*outer_len)
        o4 = Pnt( 0.5*outer_len, -0.5*outer_len, -0.5*outer_len)

    seg1 = Segment(o1,o2);
    seg2 = Segment(o2,o3);
    seg3 = Segment(o3,o4);
    seg4 = Segment(o4,o1);
    outer_face = Face(Wire([seg1,seg2,seg3,seg4]))
    
    if splitting == 8:
        outer_body = outer_face.Extrude( 0.5*outer_len*Z )
    else:
        outer_body = outer_face.Extrude( outer_len*Z )
        
    a1a = Pnt(coil_rad_outer,0,0)
    a1b = Pnt(0,coil_rad_outer,0)
    a1c = Pnt(coil_rad_outer*sqrt(2)/2,coil_rad_outer*sqrt(2)/2,0)

    b1a = Pnt(coil_rad_inner,0,0)
    b1b = Pnt(0,coil_rad_inner,0)
    b1c = Pnt(coil_rad_inner*sqrt(2)/2,coil_rad_inner*sqrt(2)/2,0)

    seg1 = ArcOfCircle(a1a, a1c, a1b)
    seg2 = Segment(a1b,b1b)
    seg3 = ArcOfCircle(b1b, b1c, b1a)
    seg4 = Segment(b1a,a1a)
    coil_face = Face(Wire([seg1,seg2,seg3,seg4]))

    coil_body1 = coil_face.Extrude( 0.5*coil_z_len*Z )
    coil_body2 = coil_body1.Move((0,0,-0.5*coil_z_len))
    coil_body3 = coil_body1.Rotate( Axis((0,0,0.25*coil_z_len), Z), 90)
    coil_body4 = coil_body2.Rotate( Axis((0,0,0.25*coil_z_len), Z), 90)
    coil_body5 = coil_body1.Rotate( Axis((0,0,0.25*coil_z_len), Z), 180)
    coil_body6 = coil_body2.Rotate( Axis((0,0,0.25*coil_z_len), Z), 180)
    coil_body7 = coil_body1.Rotate( Axis((0,0,0.25*coil_z_len), Z), 270)
    coil_body8 = coil_body2.Rotate( Axis((0,0,0.25*coil_z_len), Z), 270)

    if splitting == 8:
        coil_body = coil_body1
    elif splitting == 4:
        coil_body = coil_body1 + coil_body2
    elif splitting == 2:
        coil_body = coil_body1 + coil_body2 + coil_body3 + coil_body4
    else:
        coil_body = coil_body1 + coil_body2 + coil_body3 + coil_body4 + coil_body5 + coil_body6 + coil_body7 + coil_body8

    c1 = Pnt(0,0,core_air_z+2*core_d)
    c2 = Pnt(0,0,0)
    c3 = Pnt(2*core_d,0,0)
    c4 = Pnt(2*core_d,0,core_air_z)
    c5 = Pnt(2*core_d+core_air_x,0,core_air_z)
    c6 = Pnt(2*core_d+core_air_x,0,0)
    c7 = Pnt(4*core_d+core_air_x,0,0)
    c8 = Pnt(4*core_d+core_air_x,0,core_air_z+2*core_d)
    seg1 = Segment(c1,c2);
    seg2 = Segment(c2,c3);
    seg3 = Segment(c3,c4);
    seg4 = Segment(c4,c5);
    seg5 = Segment(c5,c6);
    seg6 = Segment(c6,c7);
    seg7 = Segment(c7,c8);
    seg8 = Segment(c8,c1);
    core_face = Face(Wire([seg1,seg2,seg3,seg4,seg5,seg6,seg7,seg8]))

    core_body1 = core_face.Extrude(core_d*Y )
    core_body2 = core_body1.Mirror(Axis((0,0.5*core_d,0), X))
    core_body3 = core_body1.Mirror(Axis((0,0.5*core_d,0), Y))
    core_body4 = core_body1.Mirror(Axis((0,0.5*core_d,0), Z))
    
    core_body5 = core_body1.Move((0,-core_d,0))
    core_body6 = core_body2.Move((0,-core_d,0))
    core_body7 = core_body3.Move((0,-core_d,0))
    core_body8 = core_body4.Move((0,-core_d,0))
    
    if splitting == 8:
        core_body = core_body1
    elif splitting == 4:
        core_body = core_body1 + core_body2
    elif splitting == 2:
        core_body = core_body1 + core_body2 + core_body3 + core_body4
    else:
        core_body = core_body1 + core_body2 + core_body3 + core_body4 + core_body5 + core_body6 + core_body7 + core_body8

    outer_body.mat("air")
    
    outer_body.faces.Min(Z).name ="dOmega_bottom"
    outer_body.faces.Max(Z).name ="dOmega_up"
    outer_body.faces.Min(X).name ="dOmega_left"
    outer_body.faces.Max(X).name ="dOmega_right"
    outer_body.faces.Min(Y).name ="dOmega_forward"
    outer_body.faces.Max(Y).name ="dOmega_back"
    
    coil_body.mat("coil")
    core_body.mat("core")
    
    if splitting == 8:
        coil_body.faces.Min(Y).name ="dOmega_forward"
        coil_body.faces.Min(X).name ="dOmega_left"
        coil_body.faces.Min(Z).name ="dOmega_bottom"
        core_body.faces.Min(Y).name ="dOmega_forward"
        core_body.faces.Min(X).name ="dOmega_left"
        core_body.faces[1].name ="dOmega_bottom"
        core_body.faces[5].name ="dOmega_bottom"
    elif splitting == 4:
        coil_body.faces.Min(Y).name ="dOmega_forward"
        coil_body.faces.Min(X).name ="dOmega_left"
        core_body.faces.Min(Y).name ="dOmega_forward"
        core_body.faces.Min(X).name ="dOmega_left"
    elif splitting == 2:
        coil_body.faces[1].name ="dOmega_forward"
        coil_body.faces[3].name ="dOmega_forward"
        core_body.faces.Min(Y).name ="dOmega_forward"

    outer_body -= coil_body
    outer_body -= core_body
    
    domains = []

    domains.append(coil_body)
    domains.append(core_body)
    domains.append(outer_body)
    
    geo = OCCGeometry(Glue(domains))
    return geo

def geometryPaperRen():
    curved = True

    outer_height = 70*10**(-3)
    coil_z_height = 20*10**(-3)

    domains = []

    o1 = Pnt(60*10**(-3),80*10**(-3),-60*10**(-3))
    o2 = Pnt(-80*10**(-3),80*10**(-3),-60*10**(-3))
    o3 = Pnt(-80*10**(-3),0,-60*10**(-3))
    o4 = Pnt(60*10**(-3),0,-60*10**(-3))
    seg1 = Segment(o1,o2);
    seg2 = Segment(o2,o3);
    seg3 = Segment(o3,o4);
    seg4 = Segment(o4,o1);
    outer_face = Face(Wire([seg1,seg2,seg3,seg4]))
    outer_body = outer_face.Extrude( outer_height*Z )

    a1a = Pnt(32.5*10**(-3),10*10**(-3),-10*10**(-3))
    a1b = Pnt(10*10**(-3),32.5*10**(-3),-10*10**(-3))
    a1c = Pnt((12.5+(20*sqrt(2)/2))*10**(-3),(12.5+(20*sqrt(2)/2))*10**(-3),-10*10**(-3))

    a2a = Pnt(-10*10**(-3),32.5*10**(-3),-10*10**(-3))
    a2b = Pnt(-32.5*10**(-3),10*10**(-3),-10*10**(-3))
    a2c = Pnt(-(12.5+(20*sqrt(2)/2))*10**(-3),(12.5+(20*sqrt(2)/2))*10**(-3),-10*10**(-3))

    a3 = Pnt(-32.5*10**(-3),0,-10*10**(-3))
    a4 = Pnt(32.5*10**(-3),0,-10*10**(-3))

    b1a = Pnt(12.5*10**(-3),10*10**(-3),-10*10**(-3))
    b1b = Pnt(10*10**(-3),12.5*10**(-3),-10*10**(-3))
    b1c = Pnt((10+(2.5*sqrt(2)/2))*10**(-3),(10+(2.5*sqrt(2)/2))*10**(-3),-10*10**(-3))

    b2a = Pnt(-12.5*10**(-3),10*10**(-3),-10*10**(-3))
    b2b = Pnt(-10*10**(-3),12.5*10**(-3),-10*10**(-3))
    b2c = Pnt((-10-(2.5*sqrt(2)/2))*10**(-3),(10+(2.5*sqrt(2)/2))*10**(-3),-10*10**(-3))

    b3 = Pnt(-12.5*10**(-3),0,-10*10**(-3))
    b4 = Pnt(12.5*10**(-3),0,-10*10**(-3))

    seg1 = Segment(b1a,a1a)
    seg2 = Segment(a1a,a4)
    seg3 = Segment(a4,b4)
    seg4 = Segment(b4,b1a)
    right_side_face = Face(Wire([seg1,seg2,seg3,seg4]))
    right_side_body = right_side_face.Extrude( coil_z_height*Z )

    if curved:
        seg1 = ArcOfCircle(a1a, a1c, a1b)
    else:
        seg1 = Segment(a1a,a1b)
    seg2 = Segment(a1b,b1b)
    if curved:
        seg3 = ArcOfCircle(b1b, b1c, b1a)
    else:
        seg3 = Segment(b1b,b1a)
    seg4 = Segment(b1a,a1a)
    up_right_corner_face = Face(Wire([seg1,seg2,seg3,seg4]))
    up_right_corner_body = up_right_corner_face.Extrude( coil_z_height*Z )

    seg1 = Segment(b1b,b2b)
    seg2 = Segment(b2b,a2a)
    seg3 = Segment(a2a,a1b)
    seg4 = Segment(a1b,b1b)
    up_side_face = Face(Wire([seg1,seg2,seg3,seg4]))
    up_side_body = up_side_face.Extrude( coil_z_height*Z )

    if curved:
        seg1 = ArcOfCircle(a2a, a2c, a2b)
    else:
        seg1 = Segment(a2a,a2b)
    seg2 = Segment(a2b,b2a)
    if curved:
        seg3 = ArcOfCircle(b2a, b2c, b2b)
    else:
        seg3 = Segment(b2a,b2b)
    seg4 = Segment(b2b,a2a)
    up_left_corner_face = Face(Wire([seg1,seg2,seg3,seg4]))
    up_left_corner_body = up_left_corner_face.Extrude( coil_z_height*Z )

    seg1 = Segment(b2a,b3)
    seg2 = Segment(b3,a3)
    seg3 = Segment(a3,a2b)
    seg4 = Segment(a2b,b2a)
    left_side_face = Face(Wire([seg1,seg2,seg3,seg4]))
    left_side_body = left_side_face.Extrude( coil_z_height*Z )

    c1 = Pnt(10*10**(-3),0,-32.5*10**(-3))
    c2 = Pnt(10*10**(-3),0,10*10**(-3))
    c3 = Pnt(-10*10**(-3),0,10*10**(-3))
    c4 = Pnt(-10*10**(-3),0,-12.5*10**(-3))
    c5 = Pnt(-35*10**(-3),0,-12.5*10**(-3))
    c6 = Pnt(-35*10**(-3),0,-1.25*10**(-3))
    c7 = Pnt(-55*10**(-3),0,-1.25*10**(-3))
    c8 = Pnt(-55*10**(-3),0,-32.5*10**(-3))
    seg1 = Segment(c1,c2);
    seg2 = Segment(c2,c3);
    seg3 = Segment(c3,c4);
    seg4 = Segment(c4,c5);
    seg5 = Segment(c5,c6);
    seg6 = Segment(c6,c7);
    seg7 = Segment(c7,c8);
    seg8 = Segment(c8,c1);
    core_face = Face(Wire([seg1,seg2,seg3,seg4,seg5,seg6,seg7,seg8]))
    core_body = core_face.Extrude( 0.5*coil_z_height*Y )

    cc1 = Pnt(-35*10**(-3),0,1.25*10**(-3))
    cc2 = Pnt(-35*10**(-3),0,10*10**(-3))
    cc3 = Pnt(-55*10**(-3),0,10*10**(-3))
    cc4 = Pnt(-55*10**(-3),0,1.25*10**(-3))
    seg1 = Segment(cc1,cc2);
    seg2 = Segment(cc2,cc3);
    seg3 = Segment(cc3,cc4);
    seg4 = Segment(cc4,cc1);
    center_core_face = Face(Wire([seg1,seg2,seg3,seg4]))
    center_core_body = center_core_face.Extrude( 0.5*coil_z_height*Y )

    outer_body.mat("outer_body")
    outer_body.faces.Min(Z).name ="bottom_face"
    outer_body.faces.Max(Z).name ="up_face"
    outer_body.faces.Min(X).name ="left_face"
    outer_body.faces.Max(X).name ="right_face"
    outer_body.faces.Min(Y).name ="forward_face"
    outer_body.faces.Max(Y).name ="back_face"

    right_side_body.mat("right_side_body")
    right_side_body.faces.Min(Y).name ="forward_face"

    up_right_corner_body.mat("up_right_corner_body")
    up_side_body.mat("up_side_body")
    up_left_corner_body.mat("up_left_corner_body")

    left_side_body.mat("left_side_body")
    left_side_body.faces.Min(Y).name ="forward_face"

    core_body.mat("core_body")
    core_body.faces.Min(Y).name ="forward_face"

    center_core_body.mat("center_core_body")
    center_core_body.faces.Min(Y).name ="forward_face"

    core_body.faces.Max(Z).name ="up_face"
    up_right_corner_body.faces.Max(Z).name ="up_face"
    up_side_body.faces.Max(Z).name ="up_face"
    up_left_corner_body.faces.Max(Z).name ="up_face"
    center_core_body.faces.Max(Z).name ="up_face"
    right_side_body.faces.Max(Z).name ="up_face"
    left_side_body.faces.Max(Z).name ="up_face"

    outer_body -= right_side_body
    outer_body -= up_right_corner_body
    outer_body -= up_side_body
    outer_body -= up_left_corner_body
    outer_body -= left_side_body
    outer_body -= core_body
    outer_body -= center_core_body

    domains.append(outer_body)
    domains.append(right_side_body)
    domains.append(up_right_corner_body)
    domains.append(up_side_body)
    domains.append(up_left_corner_body)
    domains.append(left_side_body)
    domains.append(core_body)
    domains.append(center_core_body)


    geo = OCCGeometry(Glue(domains))
    return geo



# ============================================================================
# TEAM Problem 15: Eddy current problem with conducting plate and slot
# ============================================================================

"""
TEAM Problem 15: Measurement of impedance change due to a surface-breaking slot
in a conducting plate, using an excitation coil.

Reference: TEAM Workshop Benchmark Problem 15

Two test experiments are defined:
  - Test Experiment No. 1: f = 900 Hz, N = 3790 turns
  - Test Experiment No. 2: f = 7000 Hz, N = 408 turns
"""

from netgen.occ import (
    Pnt,
    Z,
    Box,
    Sphere,
    Cylinder,
    Glue,
    OCCGeometry,
)

from netgen.webgui import Draw as DrawGeo


def team15_geometry(
    EXPERIMENT=1,
    x_coil=0.0,
    air_radius=100.0,
    symmetry=1,
    draw=False,
    h_outer=40.0,
    h_coil=10.0,
    h_plate=7.5,
    h_air = 20.0,
    h_air_inclusion=2.0,
    h_lateral_plate = 10.0
):
    """
    Build the OCC geometry of a coil above a slotted copper plate.

    All length and mesh-size parameters are specified in mm.

    Parameters
    ----------
    EXPERIMENT : int
        Benchmark configuration (1 or 2).
    x_coil : float
        Coil-axis position along x [mm].
    air_radius : float
        Radius of the spherical air domain [mm].
    symmetry : int
        1: full geometry
        2: retain y >= 0
        4: retain x >= 0 and y >= 0 (requires x_coil = 0).
    draw : bool
        If True, display the OCC geometry.
    h_outer : float
        Maximum mesh size on the outer spherical boundary [mm].
    h_coil : float
        Maximum mesh size in the coil [mm].
    h_plate : float
        Maximum mesh size in the copper plate [mm].
    h_air_inclusion : float
        Target mesh size in the air near the plate and coil [mm].
    h_plate_lateral : float
        Target mesh size on the lateral faces of the plates [mm].

    Returns
    -------
    OCCGeometry
        Geometry ready for mesh generation.
    """

    # =========================================================
    # 1. Units and experiment parameters
    # =========================================================
    mm = 1e-3

    x_coil *= mm
    air_radius *= mm

    h_outer *= mm
    h_coil *= mm
    h_plate *= mm
    h_air *= mm
    h_air_inclusion *= mm
    h_lateral_plate *= mm

    if EXPERIMENT == 1:
        a2, a1, b, lift = (
            6.15 * mm,
            12.4 * mm,
            6.15 * mm,
            0.88 * mm,
        )
        N_turns, freq = 3790, 900.0

    elif EXPERIMENT == 2:
        a2, a1, b, lift = (
            9.34 * mm,
            18.4 * mm,
            9.0 * mm,
            2.03 * mm,
        )
        N_turns, freq = 408, 7000.0

    else:
        raise ValueError("EXPERIMENT must be 1 or 2.")

    if symmetry not in (1, 2, 4):
        raise ValueError("symmetry must be 1, 2, or 4.")

    if symmetry == 4 and abs(x_coil) > 1e-12:
        raise ValueError("symmetry=4 requires x_coil=0.")

    if min(h_outer, h_coil, h_plate, h_air_inclusion) <= 0:
        raise ValueError("All mesh sizes must be positive.")

    # =========================================================
    # 2. Physical and geometrical parameters
    # =========================================================
    t_plate = 12.22 * mm

    two_c = 12.60 * mm
    w_slot = 0.28 * mm
    h_slot = 5.00 * mm

    plate_L = 120.0 * mm
    tol = 1e-6 * mm

    # =========================================================
    # 3. Spherical air domain
    # =========================================================
    air = Sphere(
        Pnt(0, 0, 0),
        air_radius,
    )
    air.faces.name = "outer"

    # =========================================================
    # 4. Copper plate with slot
    # =========================================================
    plate = Box(
        Pnt(-plate_L / 2, -plate_L / 2, -t_plate),
        Pnt(plate_L / 2, plate_L / 2, 0),
    )

    slot = Box(
        Pnt(-two_c / 2, -w_slot / 2, -h_slot),
        Pnt(two_c / 2, w_slot / 2, 0),
    )

    plate = plate - slot
    plate.solids.name = "plate"
    plate.faces[0].maxh = h_lateral_plate
    plate.faces[1].maxh = h_lateral_plate
    plate.faces[2].maxh = h_lateral_plate
    plate.faces[3].maxh = h_lateral_plate
    plate.faces[4].maxh = h_lateral_plate

    # =========================================================
    # 5. Coil: annulus with axis parallel to z
    # =========================================================
    cyl_o = Cylinder(
        Pnt(x_coil, 0, lift),
        Z,
        r=a1,
        h=b,
    )

    cyl_i = Cylinder(
        Pnt(x_coil, 0, lift),
        Z,
        r=a2,
        h=b,
    )

    coil = cyl_o - cyl_i
    coil.solids.name = "coil"

    # =========================================================
    # 6. Air excludes the plate and coil
    # =========================================================
    air = air - plate - coil
    air.solids.name = "air"

    # =========================================================
    # 7. Symmetry clipping
    # =========================================================
    if symmetry in (2, 4):

        R = air_radius + 1.0 * mm

        if symmetry == 2:
            # Keep y >= 0
            clip = Box(
                Pnt(-R, 0, -R),
                Pnt(R, R, R),
            )

        else:
            # Keep x >= 0 and y >= 0
            clip = Box(
                Pnt(0, 0, -R),
                Pnt(R, R, R),
            )

        # Name clipping planes before Boolean operations.
        for face in clip.faces:
            center = face.center

            if symmetry == 4 and abs(center[0]) < tol:
                face.name = "sym_x"

            if abs(center[1]) < tol:
                face.name = "sym_y"

        plate = plate * clip
        coil = coil * clip
        air = air * clip

        plate.solids.name = "plate"
        coil.solids.name = "coil"
        air.solids.name = "air"

    # =========================================================
    # 8. Mesh size control
    # =========================================================
    air.solids.maxh = h_air
    plate.solids.maxh = h_plate
    coil.solids.maxh = h_coil

    # =========================================================
    # 9. Outer spherical boundary mesh size
    # =========================================================
    for face in air.faces:
        if face.name == "outer":
            face.maxh = h_outer

    # =========================================================
    # 10. Name slot faces without additional refinement
    # =========================================================
    for face in plate.faces:
        center = face.center

        if (
            abs(center[0]) <= two_c / 2 + tol
            and abs(center[1]) <= w_slot / 2 + tol
            and center[2] >= -h_slot - tol
        ):
            face.name = "slot"

    # =========================================================
    # 11. Local air refinement near the coil and slot
    # =========================================================
    margin = 3.0 * h_air_inclusion

    x_min = min(-two_c / 2, x_coil - a1) - margin
    x_max = max(two_c / 2, x_coil + a1) + margin

    y_min = min(-w_slot / 2, -a1) - margin
    y_max = max(w_slot / 2, a1) + margin

    z_min = -h_slot - margin
    z_max = lift + b + margin

    # This applies a mesh size to existing air boundary faces.
    # It does not guarantee volumetric refinement of the air.
    for face in air.faces:
        center = face.center

        near_plate_or_coil = (
            x_min <= center[0] <= x_max
            and y_min <= center[1] <= y_max
            and z_min <= center[2] <= z_max
        )

        if near_plate_or_coil and face.name != "outer":
            face.maxh = h_air_inclusion

    # =========================================================
    # 12. Glue and create OCC geometry
    # =========================================================
    geo = Glue([plate, coil, air])

    if draw:
        DrawGeo(geo)

    return OCCGeometry(geo)
