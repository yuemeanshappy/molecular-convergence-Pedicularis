# Detect gene tree-species tree discordance using Quartet Scores

**Objective:** Calculate and visualize quartet-based phylogenetic concordance scores (LQ-IC, QP-IC, and EQP-IC).

**Input:** A rooted species tree and a set of rooted gene trees with matching tip labels.

**Output:** Quartet scores and reformatted tree files for visualization.

**Software and scripts:** Quartet Scores v1.0.0;`reformat_quartetscores_output.pl`;`quartet_scores_visualization.ipynb`

**Example code:**
```
# Step 1: Prepare rooted gene trees
# Remove identifiers from gene tree tip labels to match the species tree.
./data/4.tree_discordance/tree-set-alabel-filtered-rooted-singlecopy.nwk 

# Step 2: Calculate quartet scores
QuartetScores \
    -r ./data/4.tree_discordance/sptree-alabel-filtered-rooted-singlecopy.nwk \
    -e ./data/4.tree_discordance/tree-set-alabel-filtered-rooted-singlecopy.nwk \
    -o ./data/4.tree_discordance/quartet_scores_output.nwk \
    -t 10

# Step 3: Reformat quartet scores for visualization
perl reformat_quartetscores_output.pl \
    ./data/4.tree_discordance/quartet_scores_output.nwk \
    > ./data/4.tree_discordance/reformatted_quartet_scores_output.nwk

# Visualize quartet scores using Jupyter Notebook
jupyter notebook quartet_scores_visualization.ipynb


```