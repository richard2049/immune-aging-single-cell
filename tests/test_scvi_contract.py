from __future__ import annotations

import unittest

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scvi


class ScviContractTests(unittest.TestCase):
    def test_scvi_training_returns_finite_latent_embedding(self) -> None:
        """Smoke-test the Immune Aging -> scvi-tools API boundary."""

        rng = np.random.default_rng(42)

        n_cells = 64
        n_genes = 40
        n_latent = 5

        # Small integer count matrix, independent of the real aging dataset.
        counts = rng.poisson(
            lam=1.5,
            size=(n_cells, n_genes),
        ).astype(np.float32)

        obs = pd.DataFrame(
            {"batch": pd.Categorical(["batch_1"] * 32 + ["batch_2"] * 32)},
            index=[f"cell_{i}" for i in range(n_cells)],
        )

        var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

        adata = ad.AnnData(
            X=sp.csr_matrix(counts),
            obs=obs,
            var=var,
        )

        scvi.settings.seed = 42

        # Mirror the current Immune Aging scVI contract:
        # layer=None, categorical batch covariate, no continuous covariates.
        scvi.model.SCVI.setup_anndata(
            adata,
            layer=None,
            categorical_covariate_keys=["batch"],
            continuous_covariate_keys=[],
        )

        model = scvi.model.SCVI(
            adata,
            n_hidden=16,
            n_latent=n_latent,
            n_layers=1,
        )

        # This is only a technical smoke test, not meaningful model training.
        model.train(
            max_epochs=1,
            accelerator="cpu",
            devices=1,
        )

        latent = model.get_latent_representation()

        self.assertTrue(model.is_trained)
        self.assertEqual(latent.shape, (n_cells, n_latent))
        self.assertTrue(np.isfinite(latent).all())


if __name__ == "__main__":
    unittest.main()
