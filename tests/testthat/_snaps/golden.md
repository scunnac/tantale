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

# golden: tell_tales() writes the same files with the same contents

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["file", "kind", "value"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
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
          "value": ["C-terminusAAAlignment.html", "C-terminusDNAAlignment.html", "N-terminusAAAlignment.html", "N-terminusDNAAlignment.html", "TALE_CDS_all_diagnostic_regions_hmmfile.out", "allRanges.gff", "annotale/ROI_00001/TALE_DNA_parts.fasta", "annotale/ROI_00001/TALE_Protein_parts.fasta", "annotale/ROI_00001/TALE_RVDs.fasta", "annotale/ROI_00001/protocol_analyze.txt", "annotale/ROI_00001/putativeTalOrf.fasta", "annotale/ROI_00002/TALE_DNA_parts.fasta", "annotale/ROI_00002/TALE_Protein_parts.fasta", "annotale/ROI_00002/TALE_RVDs.fasta", "annotale/ROI_00002/protocol_analyze.txt", "annotale/ROI_00002/putativeTalOrf.fasta", "annotale/ROI_00003/TALE_DNA_parts.fasta", "annotale/ROI_00003/TALE_Protein_parts.fasta", "annotale/ROI_00003/TALE_RVDs.fasta", "annotale/ROI_00003/protocol_analyze.txt", "annotale/ROI_00003/putativeTalOrf.fasta", "annotale/ROI_00004/TALE_DNA_parts.fasta", "annotale/ROI_00004/TALE_Protein_parts.fasta", "annotale/ROI_00004/TALE_RVDs.fasta", "annotale/ROI_00004/protocol_analyze.txt", "annotale/ROI_00004/putativeTalOrf.fasta", "arrayReport.tsv", "domainsReport.tsv", "hitsReport.gff", "hitsReport.tsv", "hmmerSearchOut.txt", "nhmmerHumanReadableOutputOfLastRun.txt", "pseudoTalCds.fasta", "putativeTalOrf.fasta", "rvdSequences.fas", "tell_tales.log"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "stable", "volatile", "volatile", "stable", "stable", "stable", "volatile"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["3a6ea7b0799a22a407d83292576dc702", "1e03893dba346ba22b200bbe4a528444", "31e630d42844a449dbe1ba16f8b87ab9", "70a5d38c28daa4997663cf2e1e164af5", "50a112fea7605d0311a3d19a03b8c47d", "09ad861f22d2f9df8a8872c228530bd2", "01023400366a145f6449e30d2dd1a0a1", "9a57a59db49f9e87f8d28758cb8dcbfe", "1be93fdd3eba94d2b587ba3d6076a2e0", "6", "f67b5d07473a1c28b1ae256b943c2e57", "2b3c077ffe37cd2bc74d03af2dcd5727", "265d63eb6569b334e6aff457d8e9b29f", "1e821d49e3d0358bd4845c7d015c4019", "6", "5635d755a5dc5650a3978304d41d74f4", "97cd7b6b3089aa177d2dcfaf22068765", "02ac192c50561316164ec2cc010d2d99", "fea39f0c647e02f3d42aa98fbdad579b", "6", "c2822710934920edb798a173c38c989c", "a999193e4d94c2193504952e31a6b423", "02668961774948e7a387ada5a9eb6afb", "0919279857ba783d05965ae74b4df49c", "6", "790de3881ca81405dc319b2f7b2646a6", "56a152d0d74b295e9e43beecfc60d455", "c8daa90fabd215845d43a2d541379f4a", "6bdc6d26716cfdc5cfc3d75872107cb8", "89037746c22a593709fae79690de32d1", "118", "2125", "d41d8cd98f00b204e9800998ecf8427e", "b17912eec646034b980d982b8de231f1", "c874f8e5137fe7a19e91841a22dff0f6", "34"]
        }
      ]
    }

# golden: the tables tell_tales() writes, column by column

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
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
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
          "value": ["array_id", "seqnames", "start", "end", "width", "strand", "nhmmerHitID", "query_name", "hitID", "seq", "codon_count", "frameshift_count"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "double", "double", "double", "character", "character", "character", "character", "character", "double", "double"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 2, 95, 96, 7, 2, 96, 3, 96, 60, 7, 3]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["6f846cd98eba762d7fab43c14cd69e61", "02e2eaa3fd8832a516f049676ff5e7b7", "0bca32401d99360a8379f1b62f47fbf2", "272874c50573d35c85a10e19f33f366e", "c1c4650cacf683fdf5c3040070f186c3", "df802697ca41bdac37c0791962366f0c", "d417d37a19e1719a1e44dfb1340d746f", "22ad94001768b9a03fe8595ee7bfee54", "07200c470c5d94c6c624408280634fe9", "482966f124a93c41d100d22c02c569c7", "414b7407b0232bf91ffb82a3be4d1b4b", "4a4760498746c50e9492c76c1492eff6"]
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
          "value": [1, 2, 3, 4]
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
          "value": ["array_id", "seqnames", "query_name", "codon_count"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "character", "double"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 96, 96, 96]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 2, 3, 6]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["6f846cd98eba762d7fab43c14cd69e61", "02e2eaa3fd8832a516f049676ff5e7b7", "33d76af3df451c56074a2284f6a5b465", "a0edcd8088b1447234c79b2959a00425"]
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
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
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
          "value": ["array_id", "OriginalSubjectName", "Start", "End", "Strand", "NumberOfHits", "ArraySeq", "AllDomains", "SeqOfRVD", "aberrantRepeat", "N.terminusAAlength", "C.terminusAAlength", "LongestOrfLength", "OrfCovOverArrayLength", "LongestORFSeq"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "double", "double", "character", "double", "character", "logical", "character", "logical", "double", "double", "double", "double", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 2, 3, 4, 2, 3, 4, 1, 4, 1, 2, 1, 4, 3, 4]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["7f814f99eeb07530461112190e314a50", "7e13e411b075546d9a826a45224c4114", "c2adf7b6f984b5df34a5530a6c6fc23c", "da18a7477bc142b5f22868878c5fefbf", "c5aef6f127e6362ab570a403991801bf", "bb50129abf7f05f3f1b44efa233a46a0", "2e9ad64b06795de5460120521f2a2a31", "f6691062ba8ebcf78d788eb2079ade6a", "72c257668cabc7ec8b6344e79933b4e8", "2cede2adaeedcffdd9d43aaa4de93aa3", "6092490e180aea1683bf8760fdd70ad7", "1cf32b612e28232434d2cab1af5d4f2c", "9d4d3fa19d30a48c6e54fa023c18627d", "6075a5c6d9b97778adefa958cb92ec8f", "5f3de1c7abe8505aa9e67271243db3c3"]
        }
      ]
    }

# golden: a tell_tales() run loads back as a tales object

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
          "value": ["array_id", "domain_type", "position_in_crd", "dna_seq", "source_directory", "position_in_array", "aa_seq", "rvd", "seqnames"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["character", "character", "integer", "character", "character", "integer", "character", "character", "character"]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 96, 96, 96, 96, 96, 96, 96, 96]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [4, 3, 27, 60, 4, 28, 47, 10, 2]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 0, 8, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["6f846cd98eba762d7fab43c14cd69e61", "33d76af3df451c56074a2284f6a5b465", "72cbeb11d2c9db64d7b61965b9c8273b", "498e59c8411f24736f12acb1cbc3122e", "6f846cd98eba762d7fab43c14cd69e61", "a3dc0bb240729c82f0b567a37bb9498c", "b741b91defff763194ad61b55ca7eb44", "92696f2cbcfd0c5725d6ee0bfdbf00dc", "02e2eaa3fd8832a516f049676ff5e7b7"]
        }
      ]
    }

# golden: tell_tales() with frameshift correction

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["file", "kind", "value"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
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
          "value": ["C-terminusAAAlignment.html", "C-terminusDNAAlignment.html", "CorrectionAlignmentAA/CorrectionAlignmentAA_ROI_00001.html", "CorrectionAlignmentAA/CorrectionAlignmentAA_ROI_00002.html", "CorrectionAlignmentAA/CorrectionAlignmentAA_ROI_00003.html", "CorrectionAlignmentAA/CorrectionAlignmentAA_ROI_00004.html", "CorrectionAlignmentDNA/CorrectionAlignmentDNA_ROI_00001.html", "CorrectionAlignmentDNA/CorrectionAlignmentDNA_ROI_00002.html", "CorrectionAlignmentDNA/CorrectionAlignmentDNA_ROI_00003.html", "CorrectionAlignmentDNA/CorrectionAlignmentDNA_ROI_00004.html", "N-terminusAAAlignment.html", "N-terminusDNAAlignment.html", "TALE_CDS_all_diagnostic_regions_hmmfile.out", "allRanges.gff", "annotale/ROI_00001/TALE_DNA_parts.fasta", "annotale/ROI_00001/TALE_Protein_parts.fasta", "annotale/ROI_00001/TALE_RVDs.fasta", "annotale/ROI_00001/protocol_analyze.txt", "annotale/ROI_00001/putativeTalOrf.fasta", "annotale/ROI_00002/TALE_DNA_parts.fasta", "annotale/ROI_00002/TALE_Protein_parts.fasta", "annotale/ROI_00002/TALE_RVDs.fasta", "annotale/ROI_00002/protocol_analyze.txt", "annotale/ROI_00002/putativeTalOrf.fasta", "annotale/ROI_00003/TALE_DNA_parts.fasta", "annotale/ROI_00003/TALE_Protein_parts.fasta", "annotale/ROI_00003/TALE_RVDs.fasta", "annotale/ROI_00003/protocol_analyze.txt", "annotale/ROI_00003/putativeTalOrf.fasta", "annotale/ROI_00004/TALE_DNA_parts.fasta", "annotale/ROI_00004/TALE_Protein_parts.fasta", "annotale/ROI_00004/TALE_RVDs.fasta", "annotale/ROI_00004/protocol_analyze.txt", "annotale/ROI_00004/putativeTalOrf.fasta", "arrayReport.tsv", "domainsReport.tsv", "hitsReport.gff", "hitsReport.tsv", "hmmerSearchOut.txt", "nhmmerHumanReadableOutputOfLastRun.txt", "pseudoTalCds.fasta", "putativeTalOrf.fasta", "rvdSequences.fas", "tell_tales.log"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "volatile", "stable", "stable", "stable", "stable", "stable", "volatile", "volatile", "stable", "stable", "stable", "volatile"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["3a6ea7b0799a22a407d83292576dc702", "1e03893dba346ba22b200bbe4a528444", "da648aa174491c6e7507aee54e7106f9", "3fc1c6d11b1f53da9b5589b9f75e57dc", "fd58b3f23df4f6c3adcfaf6e5014768b", "c86b8f8b822554656479b9181d7b2847", "f9e4d6d02ba383e8a39d18ad4f67e520", "175fccff135050279e8aec6d6d009dcc", "3071da7222f5b73ac857e1a4cc7b270c", "a501378ea0a1f4291a29b1776bb958c3", "31e630d42844a449dbe1ba16f8b87ab9", "70a5d38c28daa4997663cf2e1e164af5", "50a112fea7605d0311a3d19a03b8c47d", "41f3669bf6f3ae99d07754ce0d1b69e5", "01023400366a145f6449e30d2dd1a0a1", "9a57a59db49f9e87f8d28758cb8dcbfe", "1be93fdd3eba94d2b587ba3d6076a2e0", "6", "f67b5d07473a1c28b1ae256b943c2e57", "2b3c077ffe37cd2bc74d03af2dcd5727", "265d63eb6569b334e6aff457d8e9b29f", "1e821d49e3d0358bd4845c7d015c4019", "6", "5635d755a5dc5650a3978304d41d74f4", "97cd7b6b3089aa177d2dcfaf22068765", "02ac192c50561316164ec2cc010d2d99", "fea39f0c647e02f3d42aa98fbdad579b", "6", "c2822710934920edb798a173c38c989c", "a999193e4d94c2193504952e31a6b423", "02668961774948e7a387ada5a9eb6afb", "0919279857ba783d05965ae74b4df49c", "6", "790de3881ca81405dc319b2f7b2646a6", "32d69186e84c0f2331db7f55621d3cf5", "c8daa90fabd215845d43a2d541379f4a", "6bdc6d26716cfdc5cfc3d75872107cb8", "89037746c22a593709fae79690de32d1", "118", "2125", "d41d8cd98f00b204e9800998ecf8427e", "b17912eec646034b980d982b8de231f1", "c874f8e5137fe7a19e91841a22dff0f6", "34"]
        }
      ]
    }

