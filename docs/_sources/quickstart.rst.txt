Quickstart
==========

Installation
------------

Just pull the repo and cd into the jdna directory.

.. code:: bash

    pip install .

or

.. code:: bash

    pipenv install .

Examples
--------

.. code::

    from jdna import Sequence

    seq = Sequence("AGTCAGTA")
    print(seq)


.. code::

    seq = Sequence.random(300)


.. code::

    s1 = Sequence.random(300)
    s2 = s1.copy().reverse_complement()

    s1.print_alignment(s2)


.. code::

    s1 = Sequence.random(300).complement()


.. code::

    s1 = Sequence.random(550)
    s1.annotate(100, 300, "GFP")
    # ..
    s1.print()

::


                                                                ----------------GFP----------------
                                                                |<START
                                                                ----      -----------RFP-----------
    0         CCCAGGACTAGCGACTTTCCGTAACGCGACCTAACACCGGCCGTTCCTTCGAGCCAGGCAAATGTTACGTCACTTCCTTAGATTT
              GGGTCCTGATCGCTGAAAGGCATTGCGCTGGATTGTGGCCGGCAAGGAAGCTCGGTCCGTTTACAATGCAGTGAAGGAATCTAAA

              ------GFP------
              -----------------------------------------RFP-----------------------------------------
    85        TGAACAGCGCCGTACCCCGATATGATATTTAGATATATAGCAGTTACACTTGGGGTTGCTATGGACTTAGATCTGCTGTATGTTT
              ACTTGTCGCGGCATGGGGCTATACTATAAATCTATATATCGTCAATGTGAACCCCAACGATACCTGAATCTAGACGACATACAAA

              -----------------------------------------RFP-----------------------------------------
    170       TCTTACCTTCCGCATCAGGGGACAATTCGCCAGTAGAATTCAGTTTGTGCGTGAGAACATAAGATTGAATCCCACGCAGGCACAA
              AGAATGGAAGGCGTAGTCCCCTGTTAAGCGGTCATCTTAAGTCAAACACGCACTCTTGTATTCTAACTTAGGGTGCGTCCGTGTT

              ---------------------RFP----------------------
    255       GCAGGGCGGGCAGACTCTATAGGTCCTAAGACCCTGAGACTGCGTCCTCAAGATACAGGTTAACAATCCCCGTATGGAGCCGTTC
              CGTCCCGCCCGTCTGAGATATCCAGGATTCTGGGACTCTGACGCAGGAGTTCTATGTCCAATTGTTAGGGGCATACCTCGGCAAG

    340       TTAGCATGACCCGACAGGTGGGCTTGGCTCGCGTAAGTTGAGTGTTGCAGATACCTGCTGCTGCGCGGTCTAGGGGGAATCGCCG
              AATCGTACTGGGCTGTCCACCCGAACCGAGCGCATTCAACTCACAACGTCTATGGACGACGACGCGCCAGATCCCCCTTAGCGGC

    425       ATTTTGACGTAGGATCGGTAATGGGCAGTAAACCCGCAACTATTTTCAGCACCAGATGCAAGTTTCCCTAGAAAGCGTCATGGTT
              TAAAACTGCATCCTAGCCATTACCCGTCATTTGGGCGTTGATAAAAGTCGTGGTCTACGTTCAAAGGGATCTTTCGCAGTACCAA

    510       TGCAATCTCCTTAGGTCACAGCAAACATAGCAGCCCCTGT
              ACGTTAGAGGAATCCAGTGTCGTTTGTATCGTCGGGGACA


.. code::

    s1 = Sequence.random(300)
    s2 = s1.copy().reverse_complement()
    s3 = s1 + s2[:30] + s2[-10:] + s2[::-1]

.. code::

    from jdna import Reaction, Sequence

    Reaction.cyclic_assembly()
