# test suite for the formats


def main(args):
    import os, glob
    try:
        dirname = os.path.dirname(__file__)
    except NameError:
        dirname = "."
    files = glob.glob(os.path.join(dirname, "test_*.py"))
    files.sort()

    for file in files:
        name = os.path.basename(file)
        name = os.path.splitext(name)[0]
        print "###### running tests in", name, "###############"
        module = __import__("Martel.test." + name).test
        module = getattr(module, name)
        module.test()

if __name__ == "__main__":
    main([])
