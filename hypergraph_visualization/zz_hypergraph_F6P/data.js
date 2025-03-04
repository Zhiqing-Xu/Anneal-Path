var data = {
    "nodes": [
        { "id": "CCCCO", "group": 0 },
        { "id": "CCCCCCO", "group": 5 },
        { "id": "CCCC=CCO", "group": 4 },
        { "id": "CCCC=CC=O", "group": 3 },
        { "id": "CCCC(O)CC=O", "group": 2 },
        { "id": "CC=O", "group": 1 },
        { "id": "CCCC=O", "group": 1 },
        { "id": "CCCCCC=O", "group": 4 },
        { "id": "CC(=O)C(=O)O", "group": 0 },
        { "id": "CCCC(O)C(CO)C(=O)O", "group": 4 },
        { "id": "CCCC(O)C(C=O)CO", "group": 3 },
        { "id": "C=O", "group": 2 },
        { "id": "O=C(O)C(=O)CCO", "group": 1 },
        { "id": "CC(=O)CO", "group": 1 },
        { "id": "NCC(=O)O", "group": 2 },
        { "id": "NC(CO)C(=O)O", "group": 1 },
        { "id": "CC(=O)C(O)(O)C=O", "group": 1 },
        { "id": "CC(O)(C=O)C(=O)O", "group": 1 },
        { "id": "CCCC(O)C(CO)C(=O)C(=O)O", "group": 2 },
        { "id": "O=CCCO", "group": 2 },
        { "id": "CCC=CO", "group": 1 },
        { "id": "CCCC(O)C(C=O)C(=O)O", "group": 3 },
        { "id": "CC(O)C=O", "group": 2 },
        { "id": "O=C(O)O", "group": 2 },
        { "id": "O=C(O)CC(=O)C(=O)O", "group": 1 },
        { "id": "O=CCC(=O)O", "group": 2 },
        { "id": "CCCC=CC(=O)C(=O)O", "group": 3 },
        { "id": "CCCC(O)CC(=O)C(=O)O", "group": 2 },
        { "id": "CCCC=CC(=O)CO", "group": 3 },
        { "id": "CCCC(O)CC(=O)CO", "group": 2 },
        { "id": "CC(=O)C(O)(O)CC(=O)C(=O)O", "group": 1 },
        { "id": "O=CO", "group": 2 },
        { "id": "CC(O)C(=O)O", "group": 1 },
        { "id": "O=CC(=O)O", "group": 2 },
        { "id": "CCCC(O)CC(=O)C(O)(O)C=O", "group": 2 },
        { "id": "CCCC(CC=O)OP(=O)(O)O", "group": 3 },
        { "id": "CCC=CCC=O", "group": 4 },
        { "id": "CCC(C(=O)O)C(O)CC=O", "group": 3 },
        { "id": "CCC(C=O)C(=O)O", "group": 2 },
        { "id": "CCCCCC(=O)CO", "group": 4 },
        { "id": "CCCCCC(=O)C(=O)O", "group": 4 },
        { "id": "CCC=CCCO", "group": 5 },
        { "id": "CCC(C(=O)O)C(O)CCO", "group": 4 },
        { "id": "CCC(C=O)C(O)CCO", "group": 3 }
    ],
        "links": [
            ["CCCCCCO", "CCCC=CCO"],
            ["CCCC=CC=O", "CCCC=CCO"],
            ["CCCC=CC=O", "CCCC(O)CC=O"],
            ["CC=O", "CCCC=O", "CCCC(O)CC=O"],
            ["CC=O", "CC(=O)C(=O)O"],
            ["CCCC=O", "CCCCO"],
            ["CCCCCCO", "CCCCCC=O"],
            ["CCCCCC=O", "CCCC=CC=O"],
            ["CCCC(O)C(CO)C(=O)O", "CCCC=CCO"],
            ["CCCC(O)C(C=O)CO", "CCCC(O)C(CO)C(=O)O"],
            ["CCCC(O)C(C=O)CO", "C=O", "CCCC(O)CC=O"],
            ["C=O", "CC(=O)C(=O)O", "O=C(O)C(=O)CCO"],
            ["CC(=O)C(=O)O", "O=C(O)C(=O)CCO"],
            ["C=O", "CC=O", "CC(=O)CO"],
            ["CC(=O)C(=O)O", "CC(=O)CO"],
            ["C=O", "NCC(=O)O", "NC(CO)C(=O)O"],
            ["CC(=O)C(=O)O", "NC(CO)C(=O)O"],
            ["C=O", "CC(=O)C(=O)O", "CC(=O)C(O)(O)C=O"],
            ["CC(=O)C(=O)O", "CC(=O)C(O)(O)C=O"],
            ["C=O", "CC(=O)C(=O)O", "CC(O)(C=O)C(=O)O"],
            ["CC(=O)C(=O)O", "CC(O)(C=O)C(=O)O"],
            ["C=O", "CC(=O)C(=O)O", "CC(=O)CO"],
            ["CCCC(O)C(C=O)CO", "CCCC(O)CC=O"],
            ["CCCC(O)C(C=O)CO", "CCCC(O)C(CO)C(=O)C(=O)O"],
            ["O=C(O)C(=O)CCO", "CCCC=O", "CCCC(O)C(CO)C(=O)C(=O)O"],
            ["CCCC(O)C(C=O)CO", "CCCC=O", "O=CCCO"],
            ["O=C(O)C(=O)CCO", "O=CCCO"],
            ["CCCC=O", "CCC=CO"],
            ["CCCCO", "CCC=CO"],
            ["CC=O", "O=CCCO"],
            ["CCCC(O)C(CO)C(=O)O", "CCCC(O)C(C=O)C(=O)O"],
            ["O=C(O)O", "CCCC(O)CC=O", "CCCC(O)C(C=O)C(=O)O"],
            ["CC(O)(C=O)C(=O)O", "CC(O)C=O", "O=C(O)O"],
            ["O=C(O)CC(=O)C(=O)O", "CC(=O)C(=O)O", "O=C(O)O"],
            ["O=C(O)CC(=O)C(=O)O", "CC(=O)C(=O)O"],
            ["CCCC(O)CC=O", "CCCC(O)C(C=O)C(=O)O"],
            ["O=CCC(=O)O", "CCCC=O", "CCCC(O)C(C=O)C(=O)O"],
            ["O=CCC(=O)O", "O=C(O)CC(=O)C(=O)O"],
            ["O=CCC(=O)O", "CC=O"],
            ["CCCC=CC=O", "CCCC=CC(=O)C(=O)O"],
            ["CCCC(O)CC(=O)C(=O)O", "CCCC=CC(=O)C(=O)O"],
            ["CC(=O)C(=O)O", "CCCC=O", "CCCC(O)CC(=O)C(=O)O"],
            ["CCCC=CC=O", "CCCC(O)C(C=O)C(=O)O"],
            ["CCCC=CC(=O)CO", "CCCC=CC=O"],
            ["CCCC=CC(=O)CO", "CCCC(O)CC(=O)CO"],
            ["CCCC(O)CC(=O)CO", "CCCC=O", "CC(=O)CO"],
            ["C=O", "CCCC(O)CC(=O)CO", "CCCC(O)CC=O"],
            ["O=C(O)CC(=O)C(=O)O", "CC=O", "CC(=O)C(O)(O)CC(=O)C(=O)O"],
            ["CC(=O)C(=O)O", "CC(=O)C(O)(O)CC(=O)C(=O)O"],
            ["CC=O", "CC(O)C(=O)O", "O=CO"],
            ["CC(=O)C(=O)O", "CC(O)C(=O)O"],
            ["CC=O", "O=CC(=O)O", "CC(=O)C(O)(O)C=O"],
            ["CCCC(O)CC=O", "CCCC(O)CC(=O)C(=O)O"],
            ["CCCC(O)CC(=O)C(O)(O)C=O", "O=CC(=O)O", "CCCC(O)CC=O"],
            ["CCCC(O)CC(=O)C(O)(O)C=O", "CCCC=O", "CC(=O)C(O)(O)C=O"],
            ["CCCC=CC=O", "CCCC(CC=O)OP(=O)(O)O"],
            ["CCCC(CC=O)OP(=O)(O)O", "CCCC(O)CC=O"],
            ["CCCCCC=O", "CCC=CCC=O"],
            ["CCC(C(=O)O)C(O)CC=O", "CCC=CCC=O"],
            ["CCC(C=O)C(=O)O", "CC=O", "CCC(C(=O)O)C(O)CC=O"],
            ["CCC(C=O)C(=O)O", "CCCC=O"],
            ["CCCCCC=O", "CCCCCC(=O)CO"],
            ["CCCC=CC(=O)CO", "CCCCCC(=O)CO"],
            ["CCCCCC=O", "CCCCCC(=O)C(=O)O"],
            ["CCCCCC(=O)C(=O)O", "CCCC=CC(=O)C(=O)O"],
            ["CCCCCCO", "CCC=CCCO"],
            ["CCC=CCCO", "CCC(C(=O)O)C(O)CCO"],
            ["CCC(C=O)C(O)CCO", "CCC(C(=O)O)C(O)CCO"],
            ["CCCC=O", "CCC(C=O)C(O)CCO", "O=CCCO"],
            ["CCC(C(=O)O)C(O)CC=O", "CCC(C(=O)O)C(O)CCO"],
            ["CCC=CCCO", "CCC=CCC=O"]
        ]
}
