import os
import pytest
import json
import yaml
import pandas as pd
from src.blastn_search import (
    blastn_and_parse,
    save_results,
)


def test_save_results_should_return_as_expected():

    results = [{"test": "result"}]

    json_file = save_results(results, "json")
    yaml_file = save_results(results, "yaml")
    csv_file = save_results(results, "csv")
    tsv_file = save_results(results, "tsv")

    for file in [json_file, yaml_file, csv_file, tsv_file]:
        assert os.path.exists(file), f"{file} was not created."

    # Check if the json output file is a valid JSON
    with open(json_file, "r") as f:
        try:
            data = json.load(f)
            assert isinstance(data, list), "Output file is not a valid JSON."
        except json.JSONDecodeError:
            pytest.fail("Output file is not a valid JSON.")

    # Check if the yaml output file is a valid YAML
    with open(yaml_file, "r") as f:
        try:
            data = yaml.safe_load(f)
            assert isinstance(data, list), "Output file is not a valid YAML."
        except yaml.YAMLError:
            pytest.fail("Output file is not a valid YAML.")

    # Check if the csv and tsv output files are valid CSV and TSV respectively
    for file, sep in [(csv_file, ","), (tsv_file, "\t")]:
        try:
            df = pd.read_csv(file, sep=sep)
            assert not df.empty, f"Output file {file} is not a valid file."
        except pd.errors.ParserError:
            pytest.fail(f"Output file {file} is not a valid file.")

    # Optional: Delete the output files after test
    for file in [json_file, yaml_file, csv_file, tsv_file]:
        os.remove(file)


def test_blastn_and_parse_should_return_as_expected():

    fasta_file = "./tests/fasta_for_pytest.fas"
    hits = 5
    output_format = "json"

    output_file = blastn_and_parse(fasta_file, hits, output_format)

    assert os.path.exists(output_file), "Output file was not created."

    with open(output_file, "r") as f:
        try:
            data = json.load(f)
            assert isinstance(data, list), "Output file is not a valid JSON."
        except json.JSONDecodeError:
            pytest.fail("Output file is not a valid JSON.")

    # Check invalid file exception
    with pytest.raises(Exception):
        blastn_and_parse("nonexistent_file.fas", hits, output_format)

    # Optional: Delete the output file after test
    os.remove(output_file)
