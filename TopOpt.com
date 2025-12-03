from driverConstants import *
from driverOptimization import OptimizationAnalysis
import driverUtils, sys
options = {
    'analysisType':OPTIMIZATION,
    'applicationName':OPTIMIZATION,
    'ask_delete':OFF,
    'job':'TopOpt',
    'noGUI':None,
    'sequential':None,
    'task':'TopOpt',
    'tmpdir':'C:\\Users\\Lorenzo\\AppData\\Local\\Temp',
}
analysis = OptimizationAnalysis(options)
status = analysis.run()
sys.exit(status)
