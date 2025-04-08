######################################
# this file makes the time history analysis convergence.
# explanation:
# if the analysis converges, only one loop happens. otherwise at the time that divergence happens
# the programm enters another while to solve divergence by changing algorithm or test or by changing
# time-step only for 0.1 duration of the time. when it solved, the programm quit the second while
# and trys to solve rest of the analysis in the first while. so the parameters (including dtAnalysis
# and tolerance) changed in the second while should be set back to initial ones except for algorithmType
# and testType because the modified ones are more soutable for that specific record. the first while finishes
# in these two situations: either currentTime reaches Tmax or the analysis dosen't converge after reducing
# dtAnalysis by the amount of minStepRatio.
# the first while try to analyze the majority of the record using fast options such as big dtAnalysis,
# and when the first while couldn't converge the analysis, the second while uses more accurate options
# such as smaller dtAnalysis for the convergance purpose.
##########
# input from other files: 
# dtInput  GMfile  scaleFac Tmax logFileID
#####

# if NTHDur dose not exit it returns 0
if {[info exist NTHDur] == 0} {
	set NTHDur "off"
}

# source getMaxResp.tcl
if {$NTHDur == "on"} {
	source getDuration.tcl;
}

# # record 1-x
# set timeSeriesID 100
# set GMfile "$gmPath/$iRec-x.txt"
# timeSeries Path $timeSeriesID -dt $dtInput -filePath $GMfile -factor $scaleFac
# pattern UniformExcitation 11 1 -accel $timeSeriesID

# # record 1-y
# set timeSeriesID 200
# set GMfile "$gmPath/$iRec-y.txt"
# timeSeries Path $timeSeriesID -dt $dtInput -filePath $GMfile -factor $scaleFac
# pattern UniformExcitation 12 2 -accel $timeSeriesID

# # record 1-z
# set timeSeriesID 300
# set GMfile "$gmPath/$iRec-z.txt"
# timeSeries Path $timeSeriesID -dt $dtInput -filePath $GMfile -factor $scaleFac
# pattern UniformExcitation 13 3 -accel $timeSeriesID

# source ../generic/defineRecorders.tcl;		# since absolute acceleration recorder needs time series tag, it must be defined after time series

# define dynamic analysis parameters
wipeAnalysis

constraints Transformation;		# how it handles boundary conditions #if you use equalDOF you must use transformation
# constraints Plain;		# how it handles boundary conditions #if you use equalDOF you must use transformation
numberer RCM;					# renumber dof's to minimize band-width (optimization) # rem: Plain is only for very small model
system BandGeneral;				# how to store and solve the system of equations in the analysis # rem: bandGeneral is the most common, however for heavy models (20 story or 3D) UmPack is suggested
# system FullGeneral;				# how to store and solve the system of equations in the analysis # rem: bandGeneral is the most common, however for heavy models (20 story or 3D) UmPack is suggested
# set testType "EnergyIncr ";	# type of convergence criteria with tolerance, max iterations. if it changed during analysis we don't set it back to initial one.
set testType "NormDispIncr ";	# type of convergence criteria with tolerance, max iterations. if it changed during analysis we don't set it back to initial one.
set algorithmType "KrylovNewton"
test NormDispIncr 1.e-4 100;
algorithm $algorithmType
integrator Newmark 0.5 0.25;	# uses Newmark's average acceleration method to compute the time history
analysis Transient;				# type of analysis: transient or static

set ok 0;						# update parameter to avoid value of previous sa
set failureFlag 0;				# update parameter to avoid value of previous sa
set tolInitial 1.e-4;			# initial tolerance for a test
set dtInitial $dtInput;			# initial time step for analysis
set minStepRatio 1.e-3;		# how much the program is allowed to reduce the time-step before announcing a non-converged analysis as collapse
if {$NTHDur == "on"} {
	set startTime [clock seconds]
}
set currentTime [getTime];
set diffTime [expr $Tmax - $currentTime];	# difference between current time and maximum time (maximum time is usually duration of the record)
set nw1 0;
while {$diffTime > $dtInput && $failureFlag == 0 } {
	incr nw1
	set nw2 0
	set dtAnalysis $dtInitial;							# update the parameter after previous iteration for the rest of analysis
	set currentTime [getTime];							# update the parameter after previous iteration for the rest of analysis
	set diffTime [expr $Tmax - $currentTime];			# update diffTime to check the next conditions in while
	set numSteps [expr int($diffTime/$dtAnalysis)+1];	# number of steps for the rest of analysis
	set tol $tolInitial;								# update the parameter after previous iteration for the rest of analysis
	test NormDispIncr 1.e-4 100;							# update the parameter after previous iteration for the rest of analysis (tol might have been changed in previous iteration)
	# puts $logFileID "------------- Running: algorithmType:$algorithmType, testType:$testType, currentTime=$currentTime, diffTime=$diffTime, dtAnalysis=$dtAnalysis, tol= $tol"
	set ok [analyze $numSteps $dtAnalysis];				# ok = 0 if analysis was completed	
	set currentTime [getTime];							# update the parameter after analysis to check the conditions in the next while
	set diffTime [expr $Tmax - $currentTime];			# update the parameter after analysis to check the conditions in the next while
	
	while {$ok != 0 && $diffTime > $dtInitial} {
		
		
		# # it works oonly by Dr. Jalali OpenSEES
		# if {[getMaxResp recTags] > 0.1} {
			# set failureFlag 1
			# puts $logFileID "\n!!!!!!!!!!! max drift reached !!!!!!!!!!"
			# break
		# }

		# it works by either original OpenSEES or Dr. Jalali OpenSEES
		# source getMaxDrift.tcl;			# maximum drift of current time is calculated using this file.
		# if {$maxDisp > 0.07} {
			# set failureFlag 1
			# puts $logFileID "\n!!!!!!!!!!! max drift reachd: maxDisp = $maxDisp !!!!!!!!!!"
			# break
		# }
		
		incr nw2
		# puts $logFileID "========== nw1 = $nw1, nw2 = $nw2 =========="
		set currentTime [getTime]
		# puts $logFileID "Analysis failed at time: $currentTime. trying alternative algorithm"
		foreach algorithmType {NewtonLineSearch ModifiedNewton Newton BFGS Broyden} param {0.65 "" "" "" ""} {
			algorithm $algorithmType $param
			foreach testType {NormDispIncr} {
				test $testType $tol 100
				set currentTime [getTime];					# update time after previous iteration in foreach
				set numSteps [expr int(0.1/$dtAnalysis)]
				# puts $logFileID "trying: algorithmType:$algorithmType, testType:$testType, currentTime=$currentTime, diffTime=$diffTime, dtAnalysis=$dtAnalysis, tol= $tol"
				set ok [analyze $numSteps $dtAnalysis]
				if {$ok == 0} {
					break
				}
			}
			if {$ok == 0} {
				break
			}
		}
		if {$ok == 0} {
			set currentTime [getTime];				 	# update time after convergance to check the condition of the first while
			set diffTime [expr $Tmax - $currentTime];	# update time after convergance to check the condition of the first while			
			break
		}
		if {$ok != 0} {
			set dtAnalysis [expr $dtAnalysis/2.]
			# set tol [expr $tol*3.]
		}
		if {[expr $dtAnalysis/$dtInitial] < $minStepRatio} {
			# puts $logFileID "minStepRatio reached: the building collapsed!"
			set failureFlag 1
			break
		}
	}
}
if {$NTHDur == "on"} {
	set endTime [clock seconds]
	set SaDuration [getDuration $startTime $endTime];		# uncomment for duration
}

if {$ok != 0 } {
	# puts $logFileID "\n-------------------Analysis Interrupted-------------------"
	puts "\n-------------------Analysis Interrupted-------------------"
} else {
	# puts $logFileID "\n-------------------Analysis Completed-------------------"
	puts "\n-------------------Analysis Completed-------------------"
}
set currentTime [getTime]
# puts $logFileID "current time: [getTime]"
if {$NTHDur == "on"} {
	# puts $logFileID "\nStart time: [clock format $startTime  -format %H:%M:%S] ([clock format $startTime  -format %D])"
	# puts $logFileID "End time:   [clock format $endTime  -format %H:%M:%S] ([clock format $endTime  -format %D])"
	# puts $logFileID "Duration(min) : [expr $SaDuration/60.]";				# uncomment for duration
	puts "Duration(min) : [expr $SaDuration/60.]";				# uncomment for duration
}
# puts $logFileID "-------------------------------------------------------- \n\n\n\n"
# flush $logFileID
remove loadPattern 11
remove loadPattern 12
remove loadPattern 13
