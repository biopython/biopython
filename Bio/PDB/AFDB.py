# Bio/PDB/AFDB.py

"""
Module to access AlphaFold 3D structure predictions
"""

import os
import requests

BASE_URL = "https://alphafold.ebi.ac.uk/files"

class AFDBList:
    """
    Interface to download structures from the AlphaFold Database.
    """
    def __init__(self, out_dir="AFDB_cache"):
        """
        Initializes the AFDBList object.

        Args:
            out_dir (str): Directory to save the structure files.
        """
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)
        print(f"[INFO] Using cache directory: {self.out_dir}")

    def retrieve(self, uniprot_id: str, fmt: str = "pdb") -> str:
        """
        Download AlphaFold prediction for a given UniProt ID.

        Args:
            uniprot_id (str): UniProt ID of the protein.
            fmt (str): Either 'pdb' or 'cif'. Default is 'pdb'.

        Returns:
            str: The path to the downloaded or cached file.

        Raises:
            ValueError: If the format is not 'pdb' or 'cif'.
            FileNotFoundError: If the structure is not found on the server.
        """
        fmt = fmt.lower()
        if fmt not in ["pdb", "cif"]:
            raise ValueError("Format must be 'pdb' or 'cif'")

        ext = "pdb" if fmt == "pdb" else "cif"
        filename = f"AF-{uniprot_id}-F1-model_v4.{ext}"
        url = f"{BASE_URL}/{filename}"
        filepath = os.path.join(self.out_dir, filename)

        if os.path.exists(filepath):
            print(f"[INFO] File already cached: {filepath}")
            return filepath

        print(f"[INFO] Downloading: {url}")
        r = requests.get(url)
        if r.status_code == 200:
            with open(filepath, "wb") as f:
                f.write(r.content)
            print(f"[INFO] Saved to: {filepath}")
            return filepath
        else:
            raise FileNotFoundError(f"Structure not found for UniProt ID: {uniprot_id}")

# This part will be used for command-line execution later
if __name__ == "__main__":
    import argparse

    # 1. Set up the parser
    parser = argparse.ArgumentParser(
        description="Download AlphaFold structure predictions from the EBI."
    )

    # 2. Define the command-line arguments
    parser.add_argument(
        "--uniprot",
        type=str,
        required=True,
        help="UniProt ID(s) to download (e.g., 'P38398' or 'P38398,Q9Y6F2')."
    )
    parser.add_argument(
        "--fmt",
        type=str,
        default="pdb",
        choices=["pdb", "cif"],
        help="Output file format: 'pdb' or 'cif'. Default is 'pdb'."
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="AFDB_cache",
        help="Directory to save the structures. Default is 'AFDB_cache'."
    )

    # 3. Parse the arguments from the command line
    args = parser.parse_args()

    # 4. Use the arguments to run your logic
    afdb = AFDBList(out_dir=args.out_dir)

    # Handle multiple UniProt IDs
    uniprot_ids = [uid.strip() for uid in args.uniprot.split(',')]

    for uniprot_id in uniprot_ids:
        try:
            afdb.retrieve(uniprot_id, fmt=args.fmt)
        except (ValueError, FileNotFoundError) as e:
            print(f"[ERROR] {e}")