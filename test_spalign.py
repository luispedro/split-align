import six
import main
sam_input = '''\
@SQ	SN:OM-RGC.v1.014920257	LN:733
ERR594391.1	77	*	0	0	*	*	0	0	CTTATCTTTGGTATGTACCTTGGATTAGAAAAAGCACTGATGAGTTTGATACTCCAGACATGAGTGAAGAGGCACCTTTTTAATGGATAGAAATATTATCT	@@@DFFFFHFHF>GDHGGJJJEIFHHIIBFFG;EHGIIJHIJ9FDGIJFGGIGEEHIBBGHGGH/=@7@C@CGHIEHFHHFD>@CCACC>;;>CCDCBC@@	AS:i:0	XS:i:0
ERR594391.1	141	*	0	0	*	*	0	0	AAGTTTGAAAAATCTTTTATGTCATAATTTAATTGTGAGATAGCACGAGACAATGCTTCATCTGCCGCTTGATTTGTAAATGCAAGATAAGCAATTCTATT	B@@DDFFEHHGHGIJJIJFHHHJHJIJHJJDHHIGHHIIHIJJGGIJGHFFGGHIJIIHIJGJIJGHEIGHEFHGGHFBDFDF>@>@CCCDDD?CCCCCDD	AS:i:0	XS:i:0
ERR594391.20	81	OM-RGC.v1.016635653	10	60	54S47M	OM-RGC.v1.016518771	544	0	*	*	NM:i:2	MD:Z:5G24G16	AS:i:37	XS:i:24	SA:Z:OM-RGC.v1.016176943,458,+,35S57M9S,13,5;
ERR594391.20	2113	OM-RGC.v1.016176943	458	13	35H57M9H	OM-RGC.v1.016518771	544	0	*	*	NM:i:5	MD:Z:6C26T5T1T4G10	AS:i:32	XS:i:25	SA:Z:OM-RGC.v1.016635653,10,-,54S47M,60,2;
ERR594391.20	161	OM-RGC.v1.016518771	544	33	40M61S	OM-RGC.v1.016635653	10	0	*	*	NM:i:1	MD:Z:25T14	AS:i:35	XS:i:28
ERR594391.250000	77	*	0	0	*	*	0	0	*	*	AS:i:0	XS:i:0
ERR594391.250000	141	*	0	0	*	*	0	0	*	*	AS:i:0	XS:i:0
'''

def test_read_samgroups():
    assert len(list(main.read_samgroup(six.StringIO(sam_input)))) == 4