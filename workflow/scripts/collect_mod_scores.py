import pandas as pd
import scanpy as sc
import pyarrow.dataset as pads

params = pd.read_json("results/finalized/optimal-gep-params/larval-w1118-testes.json")

ds = pads.dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format="arrow")

df = ds.to_table().to_pandas()

tep = df.loc[(df["module"] == 27) & (df["qval"] < params["qval"][0])]

tep_genes = tep.loc[["FBgn" in x for x in tep["X1"]]]

tep_tes = tep.loc[[("FBgn" not in x) for x in tep["X1"]]]

ral517 = sc.read_h5ad("/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/scanpy/adult-ral517-testes/celltypes.h5ad")
wt= sc.read_h5ad("/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/scanpy/adult-wt-testes/celltypes.h5ad")
larv= sc.read_h5ad("/home/mlawlor/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-scrna/results/scanpy/larval-w1118-testes/celltypes.h5ad")


##
# Gene set scores
##
sc.tl.score_genes(ral517,tep["X1"], score_name="tep")
sc.tl.score_genes(ral517,tep_genes["X1"], score_name="tep_genes")
sc.tl.score_genes(ral517,tep_tes["X1"], score_name="tep_tes")


sc.tl.score_genes(wt,tep["X1"], score_name="tep")
sc.tl.score_genes(wt,tep_genes["X1"], score_name="tep_genes")
sc.tl.score_genes(wt,tep_tes["X1"], score_name="tep_tes")


sc.tl.score_genes(larv,tep["X1"], score_name="tep")
sc.tl.score_genes(larv,tep_genes["X1"], score_name="tep_genes")
sc.tl.score_genes(larv,tep_tes["X1"], score_name="tep_tes")


ral_df = ral517.obs[["clusters","tep","tep_genes","tep_tes"]]
wt_df = wt.obs[["clusters","tep","tep_genes","tep_tes"]]
larv_df = larv.obs[["clusters","tep","tep_genes","tep_tes"]]

ral_df["dataset"] = "ral517"
wt_df["dataset"] = "wt"
larv_df["dataset"] = "larval"

res = larv_df.append(ral_df).append(wt_df)

res.to_csv("results/finalized/x-dataset-comparison/mod_scores.csv.gz")

##
# TE expression
##
ral517_te_expr = pd.melt(ral517[:,tep_tes["X1"]].to_df(),var_name="feature",value_name="expression", ignore_index=False)
wt_te_expr = pd.melt(wt[:,tep_tes["X1"]].to_df(),var_name="feature",value_name="expression", ignore_index=False)
larv_te_expr = pd.melt(larv[:,tep_tes["X1"]].to_df(),var_name="feature",value_name="expression", ignore_index=False)

ral517_te_expr["dataset"] = "ral517"
wt_te_expr["dataset"] = "wt"
larv_te_expr ["dataset"] = "larval"

res_te_expr = larv_te_expr.append(ral517_te_expr).append(wt_te_expr)

res_te_expr.to_csv("results/finalized/x-dataset-comparison/te_expression.csv.gz")
