# -*- coding: mbcs -*-
# Abaqus CAE Script: Topology Optimization Airfoil


from abaqus import *
from abaqusConstants import *
import part, sketch, material, section, assembly, step, mesh, regionToolset
import optimization, job
import csv
import os
import datetime

# =====================================================
# === USER PARAMETERS =================================
# =====================================================

file_csv = r'C:\Users\Lorenzo\Desktop\Semester project\exports\slice_2.csv'   # CSV file directory

material_name = 'PLA'
density = 1240.0        # kg/m^3            
elastic_modulus = 2.5e9 # Pa
poisson = 0.33

step_name = 'Step-1'
mesh_size = 0.005  # m

opt_name = 'Stiffness_optimization'
opt_job_name = opt_name + '-Job'

# === Results directory =======================
base_output_dir = r'C:\Users\Lorenzo\Desktop\NoSpline'   

if not os.path.exists(base_output_dir):
    os.makedirs(base_output_dir)

timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
opt_output_dir = os.path.join(base_output_dir, 'TopOpt_' + timestamp)
os.makedirs(opt_output_dir)

# Topology optimization parameters
max_design_cycle = 30
volume_ratio = 0.8

# =====================================================
# === READ CSV =========================================
# =====================================================

coordinate = []
pressure_values = []

with open(file_csv, 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)
    for row in reader:
        if len(row) < 3:
            continue
        x = float(row[0])
        y = float(row[1])
        p = float(row[2])
        coordinate.append((x, y))
        pressure_values.append(p)

if len(coordinate) == 0:
    raise ValueError('CSV empty or not readable: {}'.format(file_csv))

# =====================================================
# === CREATE 2D SKETCH WITH SEGMENTS ==================
# =====================================================

model = mdb.models['Model-1']
s = model.ConstrainedSketch(name='AirfoilSketch', sheetSize=200.0)

for i in range(len(coordinate) - 1):
    s.Line(point1=coordinate[i], point2=coordinate[i + 1])

if coordinate[0] != coordinate[-1]:
    s.Line(point1=coordinate[-1], point2=coordinate[0])

# =====================================================
# === CREATE 2D DEFORMABLE PART ========================
# =====================================================

p = model.Part(name='Airfoil', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)

# =====================================================
# === CREATE HOLES =====================================
# =====================================================

x_coords = [pt[0] for pt in coordinate]
y_coords = [pt[1] for pt in coordinate]

length = max(x_coords) - min(x_coords)
y_center = sum(y_coords) / len(y_coords)

hole_radius = length * 0.01
x1 = min(x_coords) + 0.20 * length
x2 = min(x_coords) + 0.40 * length

for (xc, yc) in [(x1, y_center), (x2, y_center)]:
    s_hole = model.ConstrainedSketch(name='__hole__', sheetSize=200.0)
    s_hole.CircleByCenterPerimeter(center=(xc, yc), point1=(xc + hole_radius, yc))
    p.Cut(sketch=s_hole)
    del model.sketches['__hole__']


# =====================================================
# === CREATE MATERIAL AND SECTION =====================
# =====================================================

mat = model.Material(name=material_name)
mat.Density(table=((density,),))
mat.Elastic(table=((elastic_modulus, poisson),))

model.HomogeneousSolidSection(
    name='Section-Airfoil',
    material=material_name,
    thickness=None
)

region = (p.faces,)
p.SectionAssignment(region=region, sectionName='Section-Airfoil')

# =====================================================
# === CREATE INSTANCE ==================================
# =====================================================

a = model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)

try:
    del a.instances['Airfoil-1']
except:
    pass

inst = a.Instance(name='Airfoil-1', part=p, dependent=OFF)

try:
    a.Set(faces=inst.faces, name='TopOptSet')
except:
    pass

# =====================================================
# === STATIC STEP ======================================
# =====================================================

model.StaticStep(
    name=step_name,
    previous='Initial',
    description='General static step.',
    nlgeom=OFF
)

# =====================================================
# === MESH ============================================
# =====================================================

a.seedPartInstance(regions=(inst,), size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)

elemType1 = mesh.ElemType(elemCode=CPS4R, elemLibrary=STANDARD)
a.setElementType(regions=(inst.faces,), elemTypes=(elemType1,))

a.generateMesh(regions=(inst,))

# =====================================================
# === BOUNDARY CONDITIONS AT HOLES ====================
# =====================================================

tol = hole_radius * 0.4
points_hole1 = []
points_hole2 = []

for edge in inst.edges:
    px, py, pz = edge.pointOn[0]
    dist1 = ((px - x1)**2 + (py - y_center)**2)**0.5
    dist2 = ((px - x2)**2 + (py - y_center)**2)**0.5

    if abs(dist1 - hole_radius) < tol:
        points_hole1.append((px, py, pz))
    elif abs(dist2 - hole_radius) < tol:
        points_hole2.append((px, py, pz))

geom_edges_hole1 = None
geom_edges_hole2 = None

try:
    if points_hole1:
        geom_edges_hole1 = inst.edges.findAt(tuple(points_hole1))
    if points_hole2:
        geom_edges_hole2 = inst.edges.findAt(tuple(points_hole2))
except:
    pass

if geom_edges_hole1:
    region_hole1 = regionToolset.Region(edges=geom_edges_hole1)
    model.DisplacementBC(
        name='BC-Hole1', createStepName='Initial',
        region=region_hole1, u1=0.0, u2=0.0, ur3=0.0
    )

if geom_edges_hole2:
    region_hole2 = regionToolset.Region(edges=geom_edges_hole2)
    model.DisplacementBC(
        name='BC-Hole2', createStepName='Initial',
        region=region_hole2, u1=0.0, u2=0.0, ur3=0.0
    )

# =====================================================
# === APPLY PRESSURE LOADS BETWEEN POINTS =============
# =====================================================

num_segments = len(coordinate) - 1 if coordinate[0] != coordinate[-1] else len(coordinate)

for i in range(num_segments):
    p1_idx = i
    p2_idx = (i + 1) % len(coordinate)

    average_pressure = (pressure_values[p1_idx] + pressure_values[p2_idx]) / 2.0

    x_mid = (coordinate[p1_idx][0] + coordinate[p2_idx][0]) / 2.0
    y_mid = (coordinate[p1_idx][1] + coordinate[p2_idx][1]) / 2.0

    if i < num_segments / 2:
        average_pressure = -average_pressure
    
    try:
        edge_found = inst.edges.findAt(((x_mid, y_mid, 0.0),))
        surf_name = 'Surf-Segment-{}'.format(i)
        region_load = a.Surface(side1Edges=edge_found, name=surf_name)
        load_name = 'Load-Segment-{}'.format(i)
        model.Pressure(
            name=load_name,
            createStepName=step_name,
            region=region_load,
            distributionType=UNIFORM,
            field='',
            magnitude=average_pressure,
            amplitude=UNSET
        )
        print("  ✔ Load applied: {} | Pressure: {:.2f} Pa".format(load_name, average_pressure))
    except Exception as e:
        print("  ✘ Error in segment {}: {}".format(i, str(e)))

# =====================================================
# === TOPOLOGY OPTIMIZATION  ==========================
# =====================================================

if 'TopOptSet' not in a.sets:
    a.Set(faces=inst.faces, name='TopOptSet')

if 'TopOptSet' in a.sets:

    try:
        topo_region = a.sets['TopOptSet']

        topo_task = model.TopologyTask(
            name='Maximize Stiffness',
            region=topo_region,
            algorithm=CONDITION_BASED_OPTIMIZATION,
            freezeBoundaryConditionRegions=ON,
            freezeLoadRegions=ON
        )

        topo_task.SingleTermDesignResponse(
            name='Strain Energy', region=MODEL,
            identifier='STRAIN_ENERGY', operation=SUM
        )

        topo_task.SingleTermDesignResponse(
            name='Volume', region=MODEL,
            identifier='VOLUME', operation=SUM
        )

        topo_task.ObjectiveFunction(
            name='Minimize Strain Energy',
            objectives=((OFF, 'Strain Energy', 1.0, 0.0, ''),)
        )

        topo_task.OptimizationConstraint(
            name='Volume Constraint',
            designResponse='Volume',
            restrictionMethod=RELATIVE_EQUAL,
            restrictionValue=volume_ratio
        )

        vector1 = ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
        region_for_demold = regionToolset.Region(faces=inst.faces)

        topo_task.TopologyDemoldControl(
            name='Restrict-1',
            region=region_for_demold,
            collisionCheckRegion=DEMOLD_REGION,
            pullDirection=vector1,
            draftAngle=0.0,
            technique=AUTO
        )

        os.chdir(opt_output_dir)

        process = mdb.OptimizationProcess(
            name=opt_name,
            model='Model-1',
            task='Maximize Stiffness',
            prototypeJob=opt_job_name,
            maxDesignCycle=max_design_cycle,
            odbMergeFrequency=2,
            dataSaveFrequency=OPT_DATASAVE_SPECIFY_CYCLE,
            saveInitial=False
        )

        process.Job(
            name=opt_job_name,
            model='Model-1',
            memory=90, memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            numCpus=1, numGPUs=0
        )

        process.submit()
        process.waitForCompletion()

        odb_path = os.path.join(opt_output_dir, opt_job_name + '.odb')

        if os.path.exists(odb_path):
            session.openOdb(name=odb_path)
            print('ODB opened:', odb_path)

        print('Optimization completed. Results in:', opt_output_dir)

    except Exception as e_opt:
        print('Topology Optimization error:', e_opt)

else:
    print('TopOptSet not found')

# =====================================================
# === FINAL MESSAGE ===================================
# =====================================================

print('\n=== COMPLETED ===')
print('CSV:', file_csv)
print('Mesh:', mesh_size)
print('Volume ratio:', volume_ratio)
print('Output:', opt_output_dir)

