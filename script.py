import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import click
import os


@click.command()
@click.option("--input_csv", prompt="Input data",
              help='Path to input data.', default=os.path.join(os.path.dirname(__file__), "genomes.csv"))
@click.option('--searched_domain', prompt='Domain', default="ACCAT",
              help='Searched domain.')
@click.option("--reference_id", prompt="Reference", default="M_AURIS_656884",
              help="Reference id.")
@click.option("--out_vis", prompt="Output histogram", default=os.path.dirname(__file__),
              help='Path to output visualisation.')
@click.option("--out_csv", prompt="Output data", default=os.path.dirname(__file__),
              help="Path to the output data.")
def make_output(input_csv, out_csv, out_vis, reference_id, searched_domain):
    """Simple script that identifies common mutations and analyzes genomes for the existence of specific substrings."""

    if not all(char in ["A", "C", "N", "T", "G"] for char in list(searched_domain)):
        raise ValueError("Domain should contain only [A, N, C, T, G] characters.")

    if not os.path.isdir(out_vis):
        raise ValueError("Path to the output visualisation is not correct.")

    if not os.path.isdir(out_csv):
        raise ValueError("Path to the output csv is not correct.")

    if not os.path.exists(input_csv):
        raise ValueError("Path to the input data is not correct.")

    if not input_csv.endswith(".csv"):
        raise ValueError("Wrong extension of the input file.")

    cleaned = clean_input(input_csv)
    output = count_duplicated(cleaned)

    if reference_id not in output["id"].unique():
        raise ValueError("There's no such id in the data.")

    output = find_mutations(output, reference_id)
    output = find_domains(output, searched_domain)

    output.to_csv(out_csv+"/out.csv", header=True, index=False, mode='w')
    click.echo('Output csv saved in %s.' % out_csv)

    draw_histogram(output, out_vis)
    click.echo('Output visualisation saved in %s.' % out_vis)


def clean_input(input_csv):
    input_data = pd.read_csv(input_csv)

    if not ("id" and "genome") in input_data.columns:
        raise ValueError("Wrong column names in the input data.")

    if not all([isinstance(value, str) for value in input_data["id"]]):
        raise ValueError("Wrong content type in \"id\" column.")

    if not all([isinstance(value, str) for value in input_data["genome"]]):
        raise ValueError("Wrong content type in \"genome\" column.")

    for index, genome in enumerate(input_data["genome"]):
        if not all(char in ["A", "C", "N", "T", "G", "-"] for char in list(genome)):
            input_data.drop([index], inplace=True)
    input_data.reset_index(drop=True, inplace=True)
    return input_data


def count_duplicated(cleaned_data):
    output_df = cleaned_data.drop(['id'], axis=1).groupby(cleaned_data.drop(['id'], axis=1).columns.tolist()
                                                          ).size().reset_index().rename(columns={0: "count"})
    col = [cleaned_data["id"][cleaned_data["genome"] == gen].iloc[0] for gen in output_df["genome"]]
    output_df["id"] = col
    output_df = output_df[["id", "genome", "count"]]
    return output_df


def find_mutations(output_df, reference_id):
    mutations_column = []
    for index, row in output_df.iterrows():
        mutation = []
        our_list = list(row["genome"])
        ref_list = list(output_df["genome"][output_df["id"] == reference_id].values[0])
        i = 0
        for nuc_our, nuc_ref in zip(our_list, ref_list):
            i += 1
            if (nuc_our in ["N", "-"]) or (nuc_ref in ["N", "-"]) or (nuc_ref == nuc_our):
                continue
            else:
                mutation.append((nuc_ref, i, nuc_our))
        mutations_column.append(mutation)
    output_df["mutations"] = mutations_column
    return output_df


def find_domains(output_df, searched_domain):
    is_there = []
    for index, row in output_df.iterrows():
        w = "".join(["(" + i + "|N)" for i in searched_domain])
        runs = re.split(r"[-]", row["genome"])
        dna = "".join(runs)

        if re.search(w, row["genome"]):
            is_there.append(1)
        elif re.search(searched_domain, dna):
            is_there.append(1)
        else:
            is_there.append(0)

    output_df["isDomainPresent"] = is_there
    return output_df


def draw_histogram(output_df, out_vis):
    mutations_column = output_df["mutations"].to_list()
    histo_data = []
    for mutations in mutations_column:
        if len(mutations) == 0:
            continue
        for mutation in mutations:
            histo_data.append(mutation[0] + str(mutation[1]) + mutation[2])

    u = np.unique(histo_data[0:10])
    ind = [np.where(u == n)[0][0] for n in histo_data[0:10]]
    bins = range(0, len(u) + 1)
    hist, bins = np.histogram(ind, bins, density=True)

    plt.bar(bins[:-1], hist, align="center", color=(0.2, 0.4, 0.6, 0.6), label="10 sequences")
    plt.xticks(np.unique(ind), [str(n) for n in set(histo_data[0:10])])
    plt.title("Histogram of mutation frequency")
    plt.legend(frameon=False)
    plt.savefig(out_vis + "/plot.pdf", bbox_inches='tight')


if __name__ == '__main__':
    make_output()
