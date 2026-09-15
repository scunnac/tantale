# golden: the tales column contract

    {
      "type": "character",
      "attributes": {},
      "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames", "dom_code"]
    }

---

    {
      "type": "character",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames", "dom_code"]
        }
      },
      "value": ["character", "character", "integer", "character", "character", "integer", "character", "character", "character", "character"]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames", "dom_code"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "integer", "character", "character", "integer", "character", "character", "character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [955, 955, 955, 955, 955, 955, 955, 955, 955, 955]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 3, 28, 389, 44, 29, 251, 20, 4, 251]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 88, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["efc8900f7e53d5e28d5d4e101eaa131a", "8182af8d20166d6fbe20a220ec417437", "9b91bd0ae10db48f84a80a60175d851d", "7ec1b5e955b40ae16abb5db5d4933f17", "d5ee7082c7cd8dc7c0705fe47f2c7e8f", "2ce65650ccfa7e06f6bd37d83239d9e7", "2484401ae770e20476d8f784f35cbaed", "e05229cbcade43667466bb8d199cd951", "e8b44ee0ee442a9abd07eb784449b83d", "237842ca696a86009cc1bda7c0e1d55b"]
        }
      ]
    }

# golden: anomalies reported for the reference fixture

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "check", "detail"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": []
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": []
        },
        {
          "type": "character",
          "attributes": {},
          "value": []
        },
        {
          "type": "character",
          "attributes": {},
          "value": []
        }
      ]
    }

# golden: the consumer requirements table

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["fn", "requirement", "columns"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["plot.tales", "plot.tales", "repeat_to_rvd_map_distalr", "tale_parts_to_rvd", "tales_align", "tales_coded_strings", "tales_compare", "tales_domain_codes", "tales_domain_codes", "tales_rvd_strings"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["all_of", "optional", "all_of", "all_of", "any_of", "all_of", "any_of", "all_of", "optional", "all_of"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["rvd, aa_seq, domain_type", "seqnames, alignment_position", "dom_code, rvd", "rvd", "rvd, dom_code", "dom_code", "aa_seq, dna_seq", "dom_code, aa_seq", "rvd", "rvd"]
        }
      ]
    }

# golden: projections of a tales onto strings and maps

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["name", "seq"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 44]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 32]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["15131f615150f0f4c8867408ae88b91c", "956c2d56456704152b551edb93c836f9"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["name", "seq"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 44]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 43]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["15131f615150f0f4c8867408ae88b91c", "29b81fed8e388680c929725726cba842"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["dom_code", "aa_seq", "rvd"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [251, 251, 251]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [251, 251, 20]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["6ac4f0003416ce44f38e420f3b5c6011", "77adce2140d992ad969182454fe71e90", "e7f345c2a026627945d4bc1e67a50858"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["repeatID", "RVD"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [251, 251]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [251, 20]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["6ac4f0003416ce44f38e420f3b5c6011", "e7f345c2a026627945d4bc1e67a50858"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["name", "seq"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 44]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 32]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["15131f615150f0f4c8867408ae88b91c", "8239718b0ee7167d6521f3e6ac83fe43"]
        }
      ]
    }

# golden: tales_compare() on four arrays

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames", "dom_code"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "integer", "character", "character", "integer", "character", "character", "character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 96, 96, 96, 96, 96, 96, 96, 96, 96]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 3, 27, 61, 4, 28, 48, 10, 1, 48]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 8, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["4dd577879c0fa2cc1dc4e8d0f0af2b06", "c134188277ff968d5ca4dcafea3bc4b7", "01ca990f7e11bb4a2a3e496cb4d78419", "4e501196aa7e7b638dbe345dbfb3e006", "dbbc90439e5d239ff31c20c5f600f781", "e2741b5905677d79fd6409739fdb79f3", "2ff02150cdcf8034efcd6115f7222d11", "c71ba09d03f3e6720f06bc90288cd9fd", "0ea7a567e121dfed99d9a2260b749655", "8f9101eb443b54d8aa8e93ca89da16f4"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["id1", "id2", "dissim"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "double"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [2304, 2304, 2304]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [48, 48, 50]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["d09c8bcce848b9007e6c395cfa38a7c2", "27c365e50e1214b948374c996a461f8b", "80f712ff1d28a4e83e2fe19d50ea020b"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["id1", "id2", "dissim", "arlem_score", "max_length"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "double", "double", "integer"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [16, 16, 16, 16, 16]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 4, 7, 7, 3]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["d338bca565b3cf6dfcca81e76f8e29dc", "d615b959dd9b0e7951486ee8d1cd6640", "77e98ba589a9c8ec30c8cba6caf73cae", "0e8968362b88525e9fc42785df3f7c1d", "48ede9167ba358a5eeb65ae6256a7e23"]
        }
      ]
    }

# golden: tales_align() on both residue layers

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames", "dom_code", "alignment_position"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "integer", "character", "character", "integer", "character", "character", "character", "character", "integer"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 3, 27, 61, 4, 28, 48, 10, 1, 48, 29]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["4dd577879c0fa2cc1dc4e8d0f0af2b06", "c134188277ff968d5ca4dcafea3bc4b7", "01ca990f7e11bb4a2a3e496cb4d78419", "4e501196aa7e7b638dbe345dbfb3e006", "dbbc90439e5d239ff31c20c5f600f781", "e2741b5905677d79fd6409739fdb79f3", "2ff02150cdcf8034efcd6115f7222d11", "c71ba09d03f3e6720f06bc90288cd9fd", "0ea7a567e121dfed99d9a2260b749655", "64b2e7463b358d87cde82009203c211b", "18131c51712d9804da786d78eb715ba0"]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames", "dom_code", "alignment_position"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "integer", "character", "character", "integer", "character", "character", "character", "character", "integer"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 3, 27, 61, 4, 28, 48, 10, 1, 48, 28]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["4dd577879c0fa2cc1dc4e8d0f0af2b06", "c134188277ff968d5ca4dcafea3bc4b7", "01ca990f7e11bb4a2a3e496cb4d78419", "4e501196aa7e7b638dbe345dbfb3e006", "dbbc90439e5d239ff31c20c5f600f781", "e2741b5905677d79fd6409739fdb79f3", "2ff02150cdcf8034efcd6115f7222d11", "c71ba09d03f3e6720f06bc90288cd9fd", "0ea7a567e121dfed99d9a2260b749655", "64b2e7463b358d87cde82009203c211b", "83ec6bac9ae821b136d892d991640012"]
        }
      ]
    }

---

    {
      "type": "integer",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["rvd", "dom_code"]
        }
      },
      "value": [29, 28]
    }

# golden: tales_group() partitions the arrays the same way

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["name", "group"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "integer"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 44]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [44, 4]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["15131f615150f0f4c8867408ae88b91c", "de3cdb2c494ef48e231902098f011e46"]
        }
      ]
    }

# golden: the data behind plot() on a tales_msa

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["column", "type", "n", "n_distinct", "n_missing", "digest"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {},
          "value": ["array_id", "position_in_array", "dom_code", "matchConsensusRepeat", "rvd", "matchConsensusRvd", "repeatClusterId", "repeatSimVsRef", "rvdSimVsRef"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "integer", "character", "logical", "character", "logical", "character", "double", "double"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [84, 84, 84, 84, 84, 84, 84, 84, 84]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [3, 28, 31, 3, 9, 2, 13, 9, 6]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0, 6, 0, 0, 0, 0, 6]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["67f4692c2c42c564b3855b6550ccdb0c", "9d08c76d0b8ea6c10421ab6807065287", "bc49e482ad0354f6bb3c7e9de8c5c754", "96c78072cca16b6b1b1d91e1869808a3", "ec227be4e3128a8c4fb7dcd602d0d9ac", "34ed212c7abbb6c010af81e0b91d89d0", "b9912fbf7dd5153953907a9be99a1f41", "d03aab87293651b09ec4663447b44669", "1d422ac000c08b1ee9f86768dd63c5fb"]
        }
      ]
    }

# golden: the consensus of the reference alignment

    {
      "type": "character",
      "attributes": {},
      "value": ["NTERM", "NN", "ND", "NN", "NI", "NK", "NN", "HD", "NN", "NG", "NG", "N*", "HD", "N*", "HD", "NI", "NN", "HD", "NG", "HD", "HD", "HD", "NG", "NN", "HD", "HD", "NG", "CTERM"]
    }

---

    {
      "type": "character",
      "attributes": {},
      "value": [null, "62", "48", "133", "127", "155", "73", "149", "157", "152", "51", "118", "95", "178", "94", "57", "64", "177", "115", "26", "26", "37", "50", "64", "94", "26", "49", null]
    }

