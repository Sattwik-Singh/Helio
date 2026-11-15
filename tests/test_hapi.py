import numpy as np

from helios.data import hapi

META = {
    "parameters": [
        {"name": "Time", "type": "isotime", "units": "UTC", "fill": None},
        {"name": "BGSE", "type": "double", "units": "nT", "fill": "-1.0E31", "size": [3]},
        {"name": "Np", "type": "double", "units": "cm^-3", "fill": "-1.0E31"},
        {"name": "Unused", "type": "double", "units": "", "fill": "-1.0E31"},
    ]
}
BODY = (
    "2018-04-19T23:45:00.000Z,1.5,-2.5,3.0,5.2\n"
    "2018-04-19T23:46:00.000Z,-1.0E31,-1.0E31,-1.0E31,5.4\n"
    "2018-04-19T23:47:00.000Z,2.5,-3.5,4.0,-1.0E31\n"
)


def test_vector_parameters_expand_and_unrequested_ones_are_skipped():
    assert hapi.column_names(META, ["BGSE", "Np"]) == ["Time", "BGSE_0", "BGSE_1", "BGSE_2", "Np"]


def test_csv_parsing_sets_time_index_and_masks_fill_values():
    df = hapi.parse_csv(BODY, META, ["BGSE", "Np"])
    assert list(df.columns) == ["BGSE_0", "BGSE_1", "BGSE_2", "Np"]
    assert str(df.index.tz) == "UTC" and len(df) == 3
    assert df.loc["2018-04-19T23:45:00Z", "BGSE_1"] == -2.5
    assert df.iloc[1][["BGSE_0", "BGSE_1", "BGSE_2"]].isna().all()
    assert np.isnan(df.iloc[2]["Np"]) and df.iloc[1]["Np"] == 5.4


def test_empty_response_gives_empty_frame():
    assert hapi.parse_csv("", META, ["Np"]).empty
