from dataclasses import replace
import math
import pytest

from rnglab12.reliability import (Observation, RUSPredictor, exact_binomial_interval,
                                  predict, train_predictor, validate)


def observation(seed, split, x=.5, remaining=100., censor=False):
    return Observation(str(seed), split, "integration", "drift", "frozen-config",
                       20., (x,), 20., None if censor else 20 + remaining,
                       20 + remaining if censor else None)


def fitted(censored_calibration=False):
    dev = [observation(f"dev-{i}", "development", i / 31, 50 + i * 3) for i in range(32)]
    cal = [observation(f"cal-{i}", "calibration", i / 31, 50 + i * 3,
                       censored_calibration) for i in range(32)]
    return train_predictor(dev, cal, application="integration", fault="drift", configuration="frozen-config")


def test_no_postissuance_features_or_future_relabeling():
    with pytest.raises(ValueError, match="leakage"):
        replace(observation("x", "development"), features_through=21)
    model = fitted()
    row = observation("holdout", "final_holdout")
    assert predict(model, row) == predict(model, replace(row, event_at=1000.))
    assert predict(model, row) == predict(model, replace(row, event_at=None, censored_at=1000.))


def test_independent_splits_and_one_issuance_per_seed():
    dev = [observation(str(i), "development") for i in range(20)]
    cal = [observation(str(i), "calibration") for i in range(20)]
    with pytest.raises(ValueError, match="Seed leakage"):
        train_predictor(dev, cal, application="integration", fault="drift", configuration="frozen-config")
    model = fitted()
    with pytest.raises(ValueError, match="leaked"):
        validate(model, [observation("dev-0", "final_holdout")])
    with pytest.raises(ValueError, match="leaked"):
        validate(model, [observation("same", "final_holdout")] * 2)
    with pytest.raises(ValueError, match="final_holdout"):
        validate(model, [observation("cal", "calibration")])


def test_frozen_model_roundtrip_and_modification_detection(tmp_path):
    model = fitted()
    path = tmp_path / "model.json"
    model.save(path)
    with pytest.raises(FileExistsError):
        model.save(path)
    assert RUSPredictor.from_dict(model.to_dict()) == model
    changed = model.to_dict()
    changed["residual_quantile"] += 1
    with pytest.raises(ValueError, match="modified"):
        RUSPredictor.from_dict(changed)
    assert predict(replace(model, residual_quantile=99), observation("h", "final_holdout")).status == "ABSTAIN"


def test_population_ood_and_already_failed_abstention():
    model = fitted()
    row = observation("h", "final_holdout")
    assert predict(model, replace(row, application="vmc")).reason == "UNSUPPORTED_POPULATION"
    assert predict(model, replace(row, fault="unknown")).reason == "UNSUPPORTED_POPULATION"
    assert predict(model, replace(row, configuration="unfrozen")).reason == "UNSUPPORTED_POPULATION"
    assert predict(model, replace(row, features=(2.,))).reason == "OUT_OF_DISTRIBUTION"
    assert predict(model, replace(row, event_at=20.)).reason == "ALREADY_FAILED"


def test_censored_calibration_not_discarded_and_conservative():
    dev = [observation(f"d{i}", "development", i / 31, 100.) for i in range(32)]
    complete = [observation(f"c{i}", "calibration", i / 31, 100.) for i in range(32)]
    censored = [replace(row, event_at=None, censored_at=30.) for row in complete]
    args = dict(application="integration", fault="drift", configuration="frozen-config")
    full_model = train_predictor(dev, complete, **args)
    censored_model = train_predictor(dev, censored, **args)
    assert censored_model.calibration_censored == 32
    assert censored_model.residual_quantile >= full_model.residual_quantile
    row = observation("h", "final_holdout")
    assert predict(censored_model, row).lower_bound <= predict(full_model, row).lower_bound
    assert predict(censored_model, row).lower_bound <= 10. + 1e-10


def test_validation_keeps_unknown_censors_and_all_denominators():
    model = fitted()
    bound = predict(model, observation("base", "final_holdout")).lower_bound
    rows = [observation("safe", "final_holdout", remaining=bound + 1),
            observation("unsafe", "final_holdout", remaining=bound / 2),
            observation("safe-censor", "final_holdout", remaining=bound + 1, censor=True),
            observation("unknown-censor", "final_holdout", remaining=bound / 2, censor=True),
            replace(observation("missing", "final_holdout"), event_at=None),
            observation("ood", "final_holdout", x=4)]
    result = validate(model, rows)
    assert result["counts"] == {"SAFE": 2, "UNSAFE": 1, "UNKNOWN": 2, "ABSTAINED": 1}
    assert result["emitted_bounds"] == 5
    assert result["censored_count"] == 2 and result["missing_outcome_count"] == 1
    assert result["unsafe_fraction_identification_bounds"] == [.2, .6]
    envelope = result["unsafe_risk_95pct_identification_and_sampling_envelope"]
    assert envelope[0] <= .2 and envelope[1] >= .6


def test_exact_intervals_and_insufficient_data():
    assert exact_binomial_interval(0, 0) is None
    assert exact_binomial_interval(0, 20)[1] == pytest.approx(1 - .025 ** (1 / 20))
    with pytest.raises(ValueError, match="Insufficient"):
        train_predictor([], [], application="integration", fault="drift", configuration="frozen-config")
    with pytest.raises(ValueError):
        exact_binomial_interval(4, 3)
