from PyPedal import *
import copy


def main():
    options1 = {
        "pedname": "Fake Pedigree 1",
        "messages": "verbose",
        "renumber": 1,
        "pedfile": "set1.ped",
        "pedformat": "asd",
        "debug_messages": True,
    }

    options2 = copy.copy(options1)
    options2["pedname"] = "Fake Pedigree 2"
    options2["pedfile"] = "set2.ped"

    set1 = pyp_newclasses.load_pedigree(options1, debug_load=True)
    print("Animals in set1.ped:")
    print(list(set1.idmap.keys()))

    set2 = pyp_newclasses.load_pedigree(options2, debug_load=True)
    print("Animals in set2.ped:")
    print(list(set2.idmap.keys()))

    print('Testing the "+" operator...')
    added = set1 + set2
    print(list(added.idmap.keys()))

    print("=" * 80)

    options3 = copy.copy(options1)
    options3["pedname"] = "Fake Pedigree 3"
    options3["pedfile"] = "set3.ped"

    set3 = pyp_newclasses.load_pedigree(options3, debug_load=True)
    print("Animals in set3.ped:")
    print(list(set3.idmap.keys()))

    print('Testing the "+" operator...')
    added2 = set1 + set3
    print(list(added2.idmap.keys()))


if __name__ == "__main__":
    main()
