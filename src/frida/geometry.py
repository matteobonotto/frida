
from argparse import ArgumentParser, Namespace
import json
from typing import Dict, Any

def parse_args() -> Namespace:
    ...

def load_geometry(path) -> Dict[Any, Any]:
    with open(path, "r") as f:
        geo = json.load(f)
    return geo


def main(args: Namespace) -> None:
    ...




if __name__ == "__main__":
    ### load the json file with the geometry and generate the mesh
    
    args = parse_args()
    ...
