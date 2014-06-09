#!/bin/bash

# This script will sync cluster / local scripts.

# BE CAREFUL!
fswatch-run "/Users/dancook/Documents/Spring_Rotation/Variant-Caller-Pipeline/scripts/" "rsync -avz --delete --exclude=\".*\" /Users/dancook/Documents/Spring_Rotation/Variant-Caller-Pipeline/scripts/  dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/scripts/"
