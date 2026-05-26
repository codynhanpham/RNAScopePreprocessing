import os, sys, json, argparse
import requests

BASE_URL = "https://api.geneontology.org/api/ontology/term/"

def getInputs():
    parser = argparse.ArgumentParser(
        description="Add Gene Ontology descriptions to gene sets downloaded from MSigDB @ https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp",
        epilog="Example: python add_go_descriptions.py -input /path/to/gene_sets.json",
        add_help=True,
        allow_abbrev=True,
    )
    parser.add_argument(
        "-input",
        metavar="Input",
        help="Path to the gene sets .json file",
        required=True,
        type=str,
        dest="input",
    )

    return parser.parse_args()


def getGoData(go_id):
    url = BASE_URL + go_id
    response = requests.get(url)
    data = response.json()
    return data


def main():
    args = getInputs()
    if not os.path.exists(args.input):
        print(f"File not found: {args.input}")
        sys.exit(1)

    with open(args.input, "r") as f:
        data = json.load(f)

    for i, gene_set in enumerate(data):
        print(f"Processing gene set {i+1}/{len(data)}: {gene_set}")
        gene_set_data = data[gene_set]
        go_id = gene_set_data["exactSource"]
        if go_id is None or go_id == "":
            go_data = {
                "label": "",
                "definition": "",
            }
        else:
            go_data = getGoData(go_id)
            label = go_data["label"] if "label" in go_data else ""
            definition = go_data["definition"] if "definition" in go_data else ""
            go_data = {
                "label": label,
                "definition": definition,
            }
        
        print(f"\tlabel: {label}")
        print(f"\tdefinition: {definition[0:80]} ...")

        gene_set_data["label"] = go_data["label"]
        gene_set_data["definition"] = go_data["definition"]

    with open(args.input, "w") as f:
        json.dump(data, f, indent=4)

    print(f"Gene set descriptions added to {args.input}")

    return None



if __name__ == "__main__":
    main()