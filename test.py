#!/usr/bin/python
import sys
from py_symmetric_matrix import *
import numpy as np
import unittest
from __init__ import *

# 4x4
DEPENDS = np.array([0.0, 1.0, 0.5, 0.8, 0.1, 0.1])
DEPENDS_FNAME = "test.depends.npy"
RANKS = np.array([0,4,5,2,3,1])
RANKS_FNAME = "test.ranks.npy"
VARLIST = ['AAA', 'HDAC1', 'HDAC3', 'QQQ']
VARLIST_FNAME = "test.varlist.txt"
ENRICHED_PAIRS = [
  ('HDAC1', 'HDAC3'),
  ('MDM2', 'S100B'),
  ('TAF1', 'TAF15'),
  ('GRB14', 'PRKCZ'),
]

ENRICHED_HASH = set(['HDAC1,HDAC3', 'MDM2,S100B', 'TAF1,TAF15', 'GRB14,PRKCZ'])

PINA_FNAME = "test.pina.txt"

class TestEnriched(unittest.TestCase):

  def setUp(self):
    np.save(DEPENDS_FNAME, DEPENDS)
    np.save(RANKS_FNAME, RANKS)
    fp = open(PINA_FNAME, 'w')
    fp.write(PINA_TEXT)
    fp.close()
    fp = open(VARLIST_FNAME, 'w')
    fp.write('\n'.join(VARLIST))
    fp.close()
    
  def test_enriched(self):
    """Confirm that enriched file is parsed and indexed correctly."""
    p = EnrichedSet(PINA_FNAME)
    for q in p.pairs:
      self.assertTrue(q in ENRICHED_PAIRS, q)
      self.assertTrue("%s,%s" % (q[0], q[1]) in ENRICHED_HASH, q)
    self.assertEqual(len(p.pairs), len(ENRICHED_PAIRS))
    self.assertTrue(p.exists('HDAC3', 'HDAC1'))
    self.assertFalse(p.exists('HDAC3', 'TAF1'))
    self.assertRaises(AssertionError, p.exists, 'HDAC3', 'XXX')
    self.assertRaises(AssertionError, p.exists, 'HDAC3', 'HDAC3')

  def test_depend(self):
    """Confirm that DependencySet works."""
    d = DependencySet(VARLIST_FNAME)
    p = EnrichedSet(PINA_FNAME)
    d.add('test', DEPENDS_FNAME, RANKS_FNAME)
    r = d.compare('test', p)
    self.assertEqual(len(r), 1)
    self.assertEqual(r[0]['name1'], 'HDAC1')
    self.assertEqual(r[0]['name2'], 'HDAC3')
    self.assertEqual(r[0]['counted'], 1)
    self.assertEqual(r[0]['confirmed'], 1)
    self.assertEqual(r[0]['test'], DEPENDS[3])


    
PINA_TEXT = '''"ID(s) interactor A"	"ID(s) interactor B"	"Alt. ID(s) interactor A"	"Alt. ID(s) interactor B"	"Alias(es) interactor A"	"Alias(es) interactor B"	"Interaction detection method(s)"	"Publication 1st author(s)"	"Publication Identifier(s)"	"Taxid interactor A"	"Taxid interactor B"	"Interaction type(s)"	"Source database(s)"	"Interaction identifier(s)"	"Confidence value(s)"	"Experimental role(s) interactor A"	"Experimental role(s) interactor B"	"Properties interactor A"	"Properties interactor B"	"HostOrganism(s)"
uniprotkb:O15379	uniprotkb:Q13547	uniprotkb:HDAC3(gene name)	uniprotkb:HDAC1(gene name)	-	-	MI:0493(in vivo)|MI:0493(in vivo)|MI:0493(in vivo)|MI:0493(in vivo)|MI:0493(in vivo)|MI:0018(two hybrid)|MI:0018(two hybrid)|MI:0492(in vitro)|MI:0493(in vivo)	-	pubmed:10944117|pubmed:11013263|pubmed:10944117|pubmed:16569215|pubmed:17255935|pubmed:17895379|pubmed:17895379|pubmed:18206970|pubmed:10490602	taxid:9606(Homo sapiens)	taxid:9606(Homo sapiens)	-|-|-|-|-|-|-|-|-	MI:0468(hprd)|MI:0468(hprd)|MI:0468(hprd)|MI:0468(hprd)|MI:0468(hprd)|MI:0465(dip)|MI:0465(dip)|MI:0468(hprd)|MI:0468(hprd)	HPRD_8700|HPRD_39359|HPRD_39351|HPRD_39514|HPRD_39848|DIP_54438|DIP_54437|HPRD_40145|HPRD_40387	-	-|-|-|-|-|-|-|-|-	-|-|-|-|-|-|-|-|-	go:GO:0006916|go:GO:0005737|go:GO:0004407|go:GO:0042826|go:GO:0000118|go:GO:0046329|go:GO:0045786|go:GO:0010553|go:GO:0051225|go:GO:0005876|go:GO:0006350|go:GO:0003714|go:GO:0017053	go:GO:0016581|go:GO:0016580|go:GO:0006916|go:GO:0006338|go:GO:0005829|go:GO:0070932|go:GO:0070933|go:GO:0004407|go:GO:0042826|go:GO:0042802|go:GO:0043922|go:GO:0045786|go:GO:0010553|go:GO:0008284|go:GO:0010552|go:GO:0010870|go:GO:0003700|go:GO:0016566|go:GO:0006350|go:GO:0016563|go:GO:0008134	-|-|-|-|-|-|-|-|-
uniprotkb:P04271	uniprotkb:Q00987	uniprotkb:S100B(gene name)	uniprotkb:MDM2(gene name)	-	-	MI:0071(molecular sieving)|MI:0071(molecular sieving)|MI:0071(molecular sieving)|MI:0071(molecular sieving)	-	pubmed:20591429|pubmed:20591429|pubmed:20591429|pubmed:20591429	taxid:9606(Homo sapiens)	taxid:9606(Homo sapiens)	MI:0407(direct interaction)|MI:0407(direct interaction)|MI:0407(direct interaction)|MI:0407(direct interaction)	MI:0471(mint)|MI:0471(mint)|MI:0471(mint)|MI:0471(mint)	MINT_46260|MINT_76938|MINT_76949|MINT_46249	-	prey|neutral component|prey|neutral component	bait|neutral component|bait|neutral component	go:GO:0050786|go:GO:0048154|go:GO:0007409|go:GO:0005509|go:GO:0048306|go:GO:0008283|go:GO:0007417|go:GO:0007611|go:GO:0005634|go:GO:0048471|go:GO:0043123|go:GO:0042803|go:GO:0001726|go:GO:0048156|go:GO:0008270	go:GO:0017163|go:GO:0005829|go:GO:0030666|go:GO:0019899|go:GO:0045184|go:GO:0042802|go:GO:0005626|go:GO:0044419|go:GO:0043518|go:GO:0071157|go:GO:0000122|go:GO:0005730|go:GO:0005654|go:GO:0002039|go:GO:0005886|go:GO:0008284|go:GO:0045931|go:GO:0032436|go:GO:0043234|go:GO:0006461|go:GO:0031648|go:GO:0034504|go:GO:0042787|go:GO:0004842|go:GO:0008270	in vitro:-1|in vitro:-1|in vitro:-1|in vitro:-1
uniprotkb:P21675	uniprotkb:Q92804	uniprotkb:TAF1(gene name)	uniprotkb:TAF15(gene name)	-	-	MI:0493(in vivo)	-	pubmed:8663456	taxid:9606(Homo sapiens)	taxid:9606(Homo sapiens)	-	MI:0468(hprd)	HPRD_5917	-	-	-	go:GO:0005524|go:GO:0000080|go:GO:0071339|go:GO:0006368|go:GO:0051123|go:GO:0017025|go:GO:0070577|go:GO:0004402|go:GO:0044419|go:GO:0002039|go:GO:0018105|go:GO:0018107|go:GO:0010552|go:GO:0032436|go:GO:0060261|go:GO:0046777|go:GO:0004674|go:GO:0000117|go:GO:0006974|go:GO:0043565|go:GO:0003713|go:GO:0005669	go:GO:0003677|go:GO:0003723|go:GO:0000166|go:GO:0005634|go:GO:0008270	-
uniprotkb:Q05513	uniprotkb:Q14449	uniprotkb:PRKCZ(gene name)	uniprotkb:GRB14(gene name)	-	-	MI:0096(pull down)|MI:0493(in vivo)	-	pubmed:12242277|pubmed:12242277	taxid:9606(Homo sapiens)	taxid:9606(Homo sapiens)	MI:0915(physical association)|-	(biogrid)|MI:0468(hprd)	BIOGRID_71824|HPRD_2150	-	prey|-	bait|-	go:GO:0005524|go:GO:0006916|go:GO:0005829|go:GO:0005768|go:GO:0043560|go:GO:0046627|go:GO:0050732|go:GO:0031333|go:GO:0018105|go:GO:0004697|go:GO:0007165|go:GO:0008270	go:GO:0000139|go:GO:0005070|go:GO:0010008|go:GO:0005792|go:GO:0005886	unspecified:32644|-'''


if __name__ == "__main__":
  unittest.main()
