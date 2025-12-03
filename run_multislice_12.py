# -*- coding: mbcs -*-
# Script Abaqus CAE: crea airfoil, vincoli, carichi e ottimizzazione topologica
# Multi-case e multi-slice
# Autore: Lorenzo + GPT-5

from abaqus import *
from abaqusConstants import *
import part, sketch, material, section, assembly, step, mesh, regionToolset
import optimization, job
import csv
import os
import datetime
import glob

# =====================================================
# === PARAMETRI UTENTE ================================
# =====================================================
base_csv_dir = r'C:\Users\Lorenzo\Desktop\Semester project\exports\Slices'
base_output_dir = r'C:\Users\Lorenzo\Desktop\TopOpt'

nome_materiale = 'Material-Airfoil'
densita = 2700.0
modulo_elastico = 70e9
poisson = 0.33

nome_step = 'Step-1'
mesh_size = 0.005  # m

max_design_cycle = 15
volume_ratio = 0.7

opt_name_base = 'Stiffness_optimization'

# =====================================================
# === CREA CARTELLA PRINCIPALE CON TIMESTAMP =========
# =====================================================
timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
top_output_dir = os.path.join(base_output_dir, 'TopOpt_' + timestamp)
if not os.path.exists(top_output_dir):
    os.makedirs(top_output_dir)

# =====================================================
# === TROVA TUTTI I CASE =============================
# =====================================================
case_dirs = sorted([d for d in os.listdir(base_csv_dir) if os.path.isdir(os.path.join(base_csv_dir, d))])

for case in case_dirs:
    case_path = os.path.join(base_csv_dir, case)
    case_output_dir = os.path.join(top_output_dir, case)
    if not os.path.exists(case_output_dir):
        os.makedirs(case_output_dir)
    
    # Trova tutti i CSV delle slice per questo case
    slice_files = sorted(glob.glob(os.path.join(case_path, '*.csv')))

    for slice_file in slice_files:
        slice_name = os.path.splitext(os.path.basename(slice_file))[0]  # es. case0_slice0
        slice_output_dir = os.path.join(case_output_dir, slice_name)
        if not os.path.exists(slice_output_dir):
            os.makedirs(slice_output_dir)
        
        print('📌 Elaborazione:', slice_name)

        # =====================================================
        # === LEGGI CSV ======================================
        # =====================================================
        coordinate = []
        pressure_values = []
        with open(slice_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)  # salta intestazione
            for row in reader:
                if len(row) < 3:
                    continue
                x = float(row[0])
                y = float(row[1])
                p = float(row[2])
                coordinate.append((x, y))
                pressure_values.append(p)
        if len(coordinate) == 0:
            raise ValueError('CSV vuoto o non leggibile: {}'.format(slice_file))

        # =====================================================
        # === CREA MODELLO, PART, AIRFOIL ===================
        # =====================================================
        model = mdb.models['Model-1']
        s = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        s.Spline(points=coordinate)
        p_part = model.Part(name='Airfoil',
                            dimensionality=TWO_D_PLANAR,
                            type=DEFORMABLE_BODY)
        p_part.BaseShell(sketch=s)
        del model.sketches['__profile__']

        # =====================================================
        # === CREA FORI ======================================
        # =====================================================
        x_coords = [pt[0] for pt in coordinate]
        y_coords = [pt[1] for pt in coordinate]
        lunghezza = max(x_coords) - min(x_coords)
        y_centrale = sum(y_coords) / len(y_coords)
        raggio_foro = lunghezza * 0.01

        x1 = min(x_coords) + 0.20 * lunghezza
        x2 = min(x_coords) + 0.40 * lunghezza

        for (xc, yc) in [(x1, y_centrale), (x2, y_centrale)]:
            s_hole = model.ConstrainedSketch(name='__hole__', sheetSize=200.0)
            s_hole.CircleByCenterPerimeter(center=(xc, yc), point1=(xc + raggio_foro, yc))
            p_part.Cut(sketch=s_hole)
            del model.sketches['__hole__']

        # =====================================================
        # === MATERIALE E SEZIONE ===========================
        # =====================================================
        mat = model.Material(name=nome_materiale)
        mat.Density(table=((densita,),))
        mat.Elastic(table=((modulo_elastico, poisson),))
        model.HomogeneousSolidSection(
            name='Section-Airfoil', material=nome_materiale, thickness=None
        )
        region = (p_part.faces,)
        p_part.SectionAssignment(region=region, sectionName='Section-Airfoil')

        # =====================================================
        # === ASSEMBLY =======================================
        # =====================================================
        a = model.rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        try:
            del a.instances['Airfoil-1']
        except Exception:
            pass
        inst = a.Instance(name='Airfoil-1', part=p_part, dependent=OFF)
        try:
            a.Set(faces=inst.faces, name='TopOptSet')
        except Exception:
            pass

        # =====================================================
        # === STEP STATICO ===================================
        # =====================================================
        model.StaticStep(
            name=nome_step,
            previous='Initial',
            description='Step statico generale per analisi quasi-statica.',
            nlgeom=OFF
        )

        # =====================================================
        # === MESH ==========================================
        # =====================================================
        a.seedPartInstance(regions=(inst,), size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
        elemType1 = mesh.ElemType(elemCode=CPS4R, elemLibrary=STANDARD)
        a.setElementType(regions=(inst.faces,), elemTypes=(elemType1,))
        a.generateMesh(regions=(inst,))

        # =====================================================
        # === BC AI FORI ====================================
        # =====================================================
        tol = raggio_foro * 0.4
        points_hole1 = []
        points_hole2 = []

        for edge in inst.edges:
            px, py, pz = edge.pointOn[0]
            dist1 = ((px - x1)**2 + (py - y_centrale)**2)**0.5
            dist2 = ((px - x2)**2 + (py - y_centrale)**2)**0.5
            if abs(dist1 - raggio_foro) < tol:
                points_hole1.append((px, py, pz))
            elif abs(dist2 - raggio_foro) < tol:
                points_hole2.append((px, py, pz))

        geom_edges_hole1 = None
        geom_edges_hole2 = None
        try:
            if points_hole1:
                geom_edges_hole1 = inst.edges.findAt(tuple(points_hole1))
            if points_hole2:
                geom_edges_hole2 = inst.edges.findAt(tuple(points_hole2))
        except Exception:
            pass

        if geom_edges_hole1:
            region_hole1 = regionToolset.Region(edges=geom_edges_hole1)
            model.DisplacementBC(name='BC-Hole1', createStepName='Initial', region=region_hole1,
                                 u1=0.0, u2=0.0, ur3=0.0)
        if geom_edges_hole2:
            region_hole2 = regionToolset.Region(edges=geom_edges_hole2)
            model.DisplacementBC(name='BC-Hole2', createStepName='Initial', region=region_hole2,
                                 u1=0.0, u2=0.0, ur3=0.0)

        # =====================================================
        # === FORZE CONCENTRATE =============================
        # =====================================================
        for i, (x_csv, y_csv) in enumerate(coordinate):
            p_val = pressure_values[i]
            min_dist = float('inf')
            closest_node = None
            for node in inst.nodes:
                nx, ny, nz = node.coordinates
                dist = ((nx - x_csv)**2 + (ny - y_csv)**2)**0.5
                if dist < min_dist:
                    min_dist = dist
                    closest_node = node
            if closest_node is None:
                continue

            region_node = regionToolset.Region(nodes=inst.nodes.sequenceFromLabels((closest_node.label,)))
            model.ConcentratedForce(
                name='Force-{}'.format(i+1),
                createStepName=nome_step,
                region=region_node,
                cf2=p_val
            )

        # =====================================================
        # === TOPOLOGY OPTIMIZATION ===========================
        # =====================================================
        if 'TopOptSet' not in a.sets:
            a.Set(faces=inst.faces, name='TopOptSet')

        if 'TopOptSet' in a.sets:
            print('🚀 Avvio configurazione Topology Optimization per', slice_name)

            try:
                topo_region = a.sets['TopOptSet']
                model.TopologyTask(name='Maximize Stiffness', region=topo_region,
                                   algorithm=CONDITION_BASED_OPTIMIZATION)

                model.optimizationTasks['Maximize Stiffness'].SingleTermDesignResponse(
                    name='Strain Energy', region=MODEL, identifier='STRAIN_ENERGY',
                    drivingRegion=None, operation=SUM, stepOptions=()
                )
                model.optimizationTasks['Maximize Stiffness'].SingleTermDesignResponse(
                    name='Volume', region=MODEL, identifier='VOLUME',
                    drivingRegion=None, operation=SUM, stepOptions=()
                )

                model.optimizationTasks['Maximize Stiffness'].ObjectiveFunction(
                    name='Minimize Strain Energy',
                    objectives=((OFF, 'Strain Energy', 1.0, 0.0, ''), )
                )

                model.optimizationTasks['Maximize Stiffness'].OptimizationConstraint(
                    name='Volume Constraint',
                    designResponse='Volume',
                    restrictionMethod=RELATIVE_EQUAL,
                    restrictionValue=volume_ratio
                )

                vector1 = ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
                region_for_demold = regionToolset.Region(faces=inst.faces)
                model.optimizationTasks['Maximize Stiffness'].TopologyDemoldControl(
                    name='Restrict-1', region=region_for_demold,
                    collisionCheckRegion=DEMOLD_REGION, pointRegion=None,
                    csys=None, pullDirection=vector1, draftAngle=0.0, technique=AUTO
                )

                # 🟢 Imposta directory di lavoro per la slice
                os.chdir(slice_output_dir)

                process = mdb.OptimizationProcess(
                    name=opt_name_base, model='Model-1', task='Maximize Stiffness',
                    description='', prototypeJob=opt_name_base + '-Job',
                    maxDesignCycle=max_design_cycle, odbMergeFrequency=2,
                    dataSaveFrequency=OPT_DATASAVE_SPECIFY_CYCLE, saveInitial=False
                )

                job_name = opt_name_base + '-Job'
                process.Job(
                    name=job_name, model='Model-1',
                    memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                    numCpus=1, numGPUs=0
                )

                process.submit()
                process.waitForCompletion()

                odb_path = os.path.join(slice_output_dir, job_name + '.odb')
                if os.path.exists(odb_path):
                    session.openOdb(name=odb_path)
                    print('📊 File ODB aperto automaticamente:', odb_path)
                else:
                    print('⚠️ File ODB non trovato in:', slice_output_dir)

                print('✅ Slice completata:', slice_name)

            except Exception as e_opt:
                print('⚠️ Errore Topology Optimization slice', slice_name, e_opt)

print('\n✅ Tutti i case e slice completati. Risultati salvati in:', top_output_dir)
