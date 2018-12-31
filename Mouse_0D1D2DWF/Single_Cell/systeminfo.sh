#!/bin/bash
if [[ "$OSTYPE" == "linux-gnu" ]]
then
	# Linux
	echo $OSTYPE
	echo -n "OS = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; lsb_release -d 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Kernel Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; uname -a 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Python Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; python --version 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Java Version = " 2>&1| tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; java -version 2>&1 |head -n 1|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo "Kepler Version = Kepler-2.5"|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo "bioKepler Version = biokepler-1.2"|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "GCC Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; gcc -v 2>&1 | tail -n 1 2>&1  |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "ICC Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; icc -v 2>&1 | tail -n 1 2>&1  |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	lscpu 2>&1 | sed 's/:/ =/g'|tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
elif [[ "$OSTYPE" == "darwin"* ]]
then
	# Mac OSX
	echo $OSTYPE
	system_profiler SPSoftwareDataType 2>&1| sed 's/:/ =/g' |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Python Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; python --version 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Java Version = " 2>&1| tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; java -version 2>&1 |head -n 1|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo "Kepler Version = Kepler-2.5"|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo "bioKepler Version = biokepler-1.2"|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "GCC Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; gcc -v 2>&1 | tail -n 1 2>&1  |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "ICC Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; icc -v 2>&1 | tail -n 1 2>&1  |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	system_profiler SPHardwareDataType 2>&1 | sed 's/:/ =/g' |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
elif [[ "$OSTYPE" == "cygwin" ]]
then
	# POSIX compatibility layer and Linux environment emulation for Windows
	echo $OSTYPE
	systeminfo 2>&1 |sed -e '/Hotfix(s)/,+22d' |sed 's/:/ =/g' |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Python Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; python --version 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "Java Version = " 2>&1| tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; java -version 2>&1 |head -n 1|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo "Kepler Version = Kepler-2.5"|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo "bioKepler Version = biokepler-1.2"|tee -a  /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "GCC Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; gcc -v 2>&1 | tail -n 1 2>&1  |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
	echo -n "ICC Version = " 2>&1 |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt; icc -v 2>&1 | tail -n 1 2>&1  |tee -a /Users/mtj/Desktop/Workflows_paper/updated_12302018/Workflow_Kepler-master/Mouse_0D1D2DWF/Single_Cell/output0DFolder/MultiScaleWFReport.txt
else
echo "No Known Operating System"
# Unknown.
fi
