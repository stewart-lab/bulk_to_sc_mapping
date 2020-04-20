import gup

db = gup.mongo_gup.MongoDb()
fcid = db.getFcidFromFcno('204')

subnums_fcno = [
        ['586','180'],
        ['587','180'],
        ['598','181'],
        ['598','182'],
        ['600','184'],
        ['600','185'],
        ['600','186'],
        ['605','188'],
        ['660','196'],
        ['760','204'],
        ['761','204'],
        ['766','205'],
        ]

for fcno in { el[1] for el in subnums_fcno}:
    db.importFlowcell(fcid = db.getFcidFromFcno(fcno)) 
nested = [db.getAllSamples(subNum=el[0]) for el in subnums_fcno]
samples = [s for sublist in nested for s in sublist[:5]]
gup.run_gup.doBclThruFfqFromSamps(db=db, fcid=fcid, Sample=samples, readOneOnly=False)
for s in samples:
    gup.bt_rsem.BtRsem(db=db, fcid=fcid, Sample = [s], ref_seq='mm10');
    '''
    try:
    except KeyError as e:
        print(f'Sample {s} failed: {e}');
    '''

ct = gup.run_gup.doCollationCalcFromSamples(
        db = db, 
        fcid=fcid,
        Sample=samples,
        ref_seq='mm10',
        readOneOnly = True)[0]

print (ct.getMetadata())

def trouble_shooot():
    ct = gup.calc_tuple.CalcTuple(node='BtRsem', db=db, Sample='061090_0093053')
    err_cts = ct.getUpstreamErrors()
    for err_ct in err_cts:
        print(err_ct)
        print(err_ct.getMetadata())

    {s: db.getSampInfoFromSamp(s)['genome'] for s in samples}
