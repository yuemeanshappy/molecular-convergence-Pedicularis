#!/bin/bash

BATCH_SIZE=100

# Get all OG folders in natural numeric order
mapfile -t OG_FOLDERS < <(
    find . -maxdepth 1 -type d -name 'OG_*' | sort -V
)

TOTAL=${#OG_FOLDERS[@]}

echo "Total OG folders found: $TOTAL"

# Get all tracked files once.
echo "Checking which OG folders are already tracked..."

declare -A TRACKED_OG

while IFS= read -r file; do
	    # file will look like:
	        # OG_123/file.txt
		    #
		        # Extract the first directory component.
			    og="${file%%/*}"

			        if [[ "$og" == OG_* ]]; then
					        TRACKED_OG["./$og"]=1
						    fi
					    done < <(git ls-files)

					    echo "Finished checking Git."

					    # Build list of remaining folders
					    REMAINING=()
					    PROCESSED_COUNT=0

					    for og in "${OG_FOLDERS[@]}"; do
						        if [[ -n "${TRACKED_OG[$og]}" ]]; then
								        ((PROCESSED_COUNT++))
									    else
										            REMAINING+=("$og")
											        fi
											done

											REMAINING_TOTAL=${#REMAINING[@]}

											echo "Already tracked OG folders: $PROCESSED_COUNT"
											echo "Remaining OG folders: $REMAINING_TOTAL"
											echo

											if [ "$REMAINING_TOTAL" -eq 0 ]; then
												    echo "All OG folders are already tracked."
												        exit 0
											fi

											echo "Starting processing of $REMAINING_TOTAL remaining folders."
											echo

											BATCH_NUM=1

											for ((i=0; i<REMAINING_TOTAL; i+=BATCH_SIZE)); do

												    END=$((i + BATCH_SIZE))

												        if (( END > REMAINING_TOTAL )); then
														        END=$REMAINING_TOTAL
															    fi

															        echo "=========================================="
																    echo "Batch $BATCH_NUM"
																        echo "Folders $((i+1)) - $END of $REMAINING_TOTAL"
																	    echo "=========================================="

																	        # Add this batch
																		    for ((j=i; j<END; j++)); do
																			            echo "Adding ${REMAINING[j]}"
																				            git add -- "${REMAINING[j]}" || {
																						                echo "ERROR: git add failed for ${REMAINING[j]}"
																					                exit 1
																							        }
																								    done

																								        # Check whether anything was staged
																									    if git diff --cached --quiet; then
																										            echo "Nothing to commit."
																											            BATCH_NUM=$((BATCH_NUM + 1))
																												            continue
																													        fi

																														    echo
																														        echo "Creating commit..."

																															    if ! git commit -m "add data"; then
																																            echo
																																	            echo "ERROR: git commit failed."
																																		            echo "Stopping."
																																			            exit 1
																																				        fi

																																					    echo
																																					        echo "Pushing batch $BATCH_NUM..."

																																						    if ! git push origin main; then
																																							            echo
																																								            echo "ERROR: git push failed for batch $BATCH_NUM."
																																									            echo
																																										            echo "The commit was created locally."
																																											            echo "Fix the push problem and run this script again."
																																												            exit 1
																																													        fi

																																														    echo
																																														        echo "Batch $BATCH_NUM successfully pushed."
																																															    echo

																																															        BATCH_NUM=$((BATCH_NUM + 1))
																																															done

																																															echo "=========================================="
																																															echo "DONE"
																																															echo "All OG folders have been processed."
																																															echo "=========================================="
