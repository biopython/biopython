#!/usr/bin/env python
"""
Examples of using Bio.PDB.AlphaFoldList to access AlphaFold predicted structures.

The AlphaFold Protein Structure Database contains AI-predicted 3D structures
for millions of proteins. This module provides a PDBList-like interface for
downloading these predictions.

Author: BioPython Contributors
License: Biopython License
"""

from Bio.PDB.AlphaFoldList import AlphaFoldList


def example_1_basic_download():
    """Example 1: Download a single structure by UniProt ID."""
    print("=" * 60)
    print("Example 1: Download Single Structure")
    print("=" * 60)

    # Initialize AlphaFoldList
    afl = AlphaFoldList(local_dir="alphafold_structures", verbose=True)

    # Download mouse ABL1 kinase (P00520) in PDB format
    structure_file = afl.retrieve_alphafold_structure(
        "P00520",  # UniProt accession
        file_format="pdb",  # or "mmCif"
    )

    print(f"\nDownloaded to: {structure_file}")


def example_2_with_confidence_data():
    """Example 2: Download structure with confidence scores and PAE."""
    print("\n" + "=" * 60)
    print("Example 2: Download with Confidence Data")
    print("=" * 60)

    afl = AlphaFoldList(local_dir="alphafold_structures", verbose=True)

    # Download with additional quality metrics
    structure_file = afl.retrieve_alphafold_structure(
        "P12345",  # Human Aspartate aminotransferase
        file_format="mmCif",
        include_pae=True,  # Predicted Aligned Error
        include_confidence=True,  # Per-residue pLDDT scores
    )

    print(f"\nStructure and confidence data downloaded to: {structure_file}")


def example_3_query_metadata():
    """Example 3: Query structure metadata before downloading."""
    print("\n" + "=" * 60)
    print("Example 3: Query Metadata First")
    print("=" * 60)

    afl = AlphaFoldList(verbose=True)

    # Get metadata for a protein
    uniprot_id = "Q9Y6K9"  # NF-kappa-B essential modulator
    metadata = afl.get_structure_metadata(uniprot_id)

    if metadata:
        entry = metadata[0]
        print(f"\nProtein: {entry['uniprotDescription']}")
        print(f"Organism: {entry['organismScientificName']}")
        print(f"Entry ID: {entry['entryId']}")
        print(f"Latest Version: {entry['latestVersion']}")
        print(f"Confidence: {entry['globalMetricValue']:.2f}")
        print(f"Sequence Length: {entry['sequenceEnd'] - entry['sequenceStart'] + 1}")

        # Now download it
        print("\nDownloading structure...")
        afl.retrieve_alphafold_structure(uniprot_id)


def example_4_specific_version():
    """Example 4: Download a specific model version."""
    print("\n" + "=" * 60)
    print("Example 4: Download Specific Version")
    print("=" * 60)

    afl = AlphaFoldList(local_dir="alphafold_structures", verbose=True)

    # AlphaFold regularly updates predictions (v1, v2, v3, ...)
    # Download version 4 specifically
    structure_file = afl.retrieve_alphafold_structure(
        "P00520", file_format="pdb", version=4  # Specific version
    )

    print(f"\nDownloaded version 4: {structure_file}")


def example_5_bulk_download():
    """Example 5: Download multiple structures in parallel."""
    print("\n" + "=" * 60)
    print("Example 5: Bulk Download Multiple Structures")
    print("=" * 60)

    afl = AlphaFoldList(local_dir="alphafold_structures", verbose=True)

    # List of kinases to download
    kinases = [
        "P00520",  # ABL1
        "P00533",  # EGFR
        "P04629",  # NTRK1
        "P06239",  # LCK
    ]

    print(f"\nDownloading {len(kinases)} kinase structures...")
    afl.download_structures_by_accessions(
        kinases, file_format="mmCif", max_num_threads=4  # Parallel downloads
    )

    print("\nAll structures downloaded!")


def example_6_list_proteomes():
    """Example 6: List available organism proteomes."""
    print("\n" + "=" * 60)
    print("Example 6: Browse Available Proteomes")
    print("=" * 60)

    afl = AlphaFoldList(verbose=True)

    # Get list of all available organism proteomes
    proteomes = afl.get_available_proteomes()

    print(f"\nFound {len(proteomes)} available proteomes\n")

    # Display model organisms
    print(f"{'Organism':<15} {'Tax ID':<10} {'UniProt ID':<15} {'Filename'}")
    print("-" * 80)

    common_organisms = ["HUMAN", "MOUSE", "ECOLI", "YEAST", "DROME"]
    for proteome in proteomes:
        if any(org in proteome["organism"] for org in common_organisms):
            print(
                f"{proteome['organism']:<15} {proteome['taxid']:<10} "
                f"{proteome['uniprot_id']:<15} {proteome['filename']}"
            )


def example_7_download_proteome():
    """Example 7: Download complete organism proteome (WARNING: Large files!)."""
    print("\n" + "=" * 60)
    print("Example 7: Download Complete Proteome")
    print("=" * 60)
    print("\nWARNING: Proteome downloads are very large (GBs)!")
    print("This example shows the code but does not execute the download.\n")

    afl = AlphaFoldList(local_dir="alphafold_proteomes", verbose=True)

    # Example: Download E. coli proteome (smaller, ~456 MB)
    # Uncomment to actually download:
    # afl.download_proteome(
    #     "ECOLI",  # Organism code
    #     extract=True,  # Automatically extract .tar archive
    # )

    print("To download, uncomment the download_proteome() call above.")
    print("\nAvailable organisms:")
    print("- HUMAN: Human (~4.8 GB)")
    print("- MOUSE: Mouse (~3.5 GB)")
    print("- ECOLI: E. coli (~456 MB)")
    print("- YEAST: Yeast (~1.0 GB)")
    print("- And 40+ more organisms")


def example_8_integrate_with_biopython():
    """Example 8: Parse downloaded structures with Bio.PDB."""
    print("\n" + "=" * 60)
    print("Example 8: Integration with Bio.PDB Parser")
    print("=" * 60)

    from Bio.PDB import PDBParser, MMCIFParser

    afl = AlphaFoldList(local_dir="alphafold_structures", verbose=True)

    # Download structure
    pdb_file = afl.retrieve_alphafold_structure("P00520", file_format="pdb")

    if pdb_file:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ABL1", pdb_file)

        # Analyze the structure
        print(f"\nStructure ID: {structure.id}")
        print(f"Number of models: {len(structure)}")

        # Get first model
        model = structure[0]
        print(f"Number of chains: {len(model)}")

        # Get first chain
        chain = list(model.get_chains())[0]
        print(f"Chain {chain.id}: {len(chain)} residues")

        # Count atoms
        atom_count = sum(1 for _ in structure.get_atoms())
        print(f"Total atoms: {atom_count}")


def example_9_directory_organization():
    """Example 9: Different directory organization options."""
    print("\n" + "=" * 60)
    print("Example 9: Directory Organization")
    print("=" * 60)

    # Option 1: Hierarchical directory tree (default)
    print("\nOption 1: Hierarchical tree (PDB-style)")
    afl1 = AlphaFoldList(local_dir="alphafold_tree", verbose=False)
    afl1.flat_tree = False  # default
    file1 = afl1.retrieve_alphafold_structure("P00520")
    print(f"Saved to: {file1}")
    # Creates: alphafold_tree/p0/AF-P00520-F1-model_v6.cif

    # Option 2: Flat directory
    print("\nOption 2: Flat directory")
    afl2 = AlphaFoldList(local_dir="alphafold_flat", verbose=False)
    afl2.flat_tree = True
    file2 = afl2.retrieve_alphafold_structure("P00520")
    print(f"Saved to: {file2}")
    # Creates: alphafold_flat/AF-P00520-F1-model_v6.cif

    # Option 3: Custom directory for each protein
    print("\nOption 3: Custom directory")
    afl3 = AlphaFoldList(verbose=False)
    file3 = afl3.retrieve_alphafold_structure("P00520", pdir="proteins/ABL1/")
    print(f"Saved to: {file3}")
    # Creates: proteins/ABL1/AF-P00520-F1-model_v6.cif


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("Bio.PDB.AlphaFoldList Usage Examples")
    print("=" * 60)

    examples = [
        example_1_basic_download,
        example_2_with_confidence_data,
        example_3_query_metadata,
        example_4_specific_version,
        example_5_bulk_download,
        example_6_list_proteomes,
        example_7_download_proteome,
        # example_8_integrate_with_biopython,  # Requires parsing
        # example_9_directory_organization,  # Creates many files
    ]

    for i, example_func in enumerate(examples, 1):
        try:
            example_func()
        except Exception as e:
            print(f"\n[ERROR in Example {i}]: {e}")

    print("\n" + "=" * 60)
    print("Examples complete!")
    print("=" * 60)


if __name__ == "__main__":
    # Run only example 1 and 3 (lightweight examples)
    print("Running lightweight examples (1, 3, 6)...\n")
    example_1_basic_download()
    example_3_query_metadata()
    example_6_list_proteomes()

    print("\n\nTo see all examples, uncomment the main() call below:")
    print("# main()")
