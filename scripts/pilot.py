from pathlib import Path
import json
from rnglab12.execution import APPLICATIONS,run_tasks,ROOT
if __name__=='__main__':
    tasks=[{'application':a,'split':'pilot','group':i} for a in APPLICATIONS for i in range(4)]
    result=run_tasks(ROOT/'runs/pilot',tasks,workers=2)
    print(json.dumps({'tasks':len(result),'per_task_seconds':{a:sum(r['elapsed_seconds'] for r in result if r['application']==a)/4 for a in APPLICATIONS}},indent=2))
