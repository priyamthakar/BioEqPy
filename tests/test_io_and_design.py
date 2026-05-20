from pathlib import Path

import pytest

from bioeqpy.core.exceptions import ValidationError
from bioeqpy.io.loaders import infer_parameters, load_dataset, load_study


DATA = Path(__file__).parent / "reference_data" / "example_2x2.csv"


def test_load_dataset_and_infer_parameters():
    data = load_dataset(DATA)

    assert infer_parameters(data) == ["AUClast", "AUCinf", "Cmax"]
    assert set(data.columns[:4]) == {"subject", "sequence", "period", "treatment"}


def test_load_study_detects_2x2_and_logs_values():
    study = load_study(DATA, "AUClast")

    assert study.design.name == "2x2"
    assert study.design.sequences == ["RT", "TR"]
    assert study.parameter_name == "AUClast"
    assert len(study.values) == 16


def test_missing_parameter_raises():
    with pytest.raises(ValidationError):
        load_study(DATA, "AUC0_24")

