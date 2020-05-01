//
//	rtf3_recon_2 (motion correction phase)
//
//	K. Oshio	3-10-2011
//
//  === plans ===
//  4 pts correction    1-21-2013
//

#import <RecKit/RecKit.h>
#import "RecImagePW.h"	// PW category
#import <sys/time.h>
#import <RecKit/timer_macros.h>

int
main(int ac, char *av[])
{
    @autoreleasepool {  // outer
	RecLoop			*sl, *ky, *kx, *ch, *pe, *rd;
	int				i;  // iteration
    NSString        *path;

	RecImage		*img, *raw, *pw, *pws;          // multi-ch
	RecImage		*img_c, *pw_c, *pwr_c, *pws_c;  // combined
    RecImage        *img_b, *pw_b;                  // re-proj of background (stationary part)
    float           mx_b, mx_c;
    RecImage        *corr, *sft;
//    RecImage        *mask;
	RecImage		*radMap, *traj;                 // warp map, k-traj
    int             mode = 0;                       // 0: shift, 1: shift4
    RecGridder      *grid;

TIMER_ST
	img = [RecImage imageFromFile:@"img" relativePath:YES]; // input 1
//    [img pcorr];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [RecLoop findLoop:@"Slice"];
	ky = [RecLoop findLoop:@"ky"];
	kx = [RecLoop findLoop:@"kx"];
TIMER_END("read img");

    @autoreleasepool {  // inner 1
        [img pcorr];
    //    img_c = [img complexCombineForLoop:ch withCoil:GE_Card_8_new];
        img_c = [img combineForLoop:ch withCoil:GE_Card_8_new];
        [img_c saveAsKOImage:@"img_0.img"];
    // read mask (created on OsiriX)
//        mask = [RecImage imageOfType:RECIMAGE_REAL withLoops:sl, ky, kx, nil];
//        [mask initWithRawImage:@"mask.raw"];        // input 2
//    [mask saveAsKOImage:@"mask.img"];
        pw = [RecImage imageFromFile:@"pw" relativePath:YES];	// multi-ch, complex, input 3
        pe = [RecLoop findLoop:@"Phase"];
        rd = [RecLoop findLoop:@"rd_zf"];
    TIMER_END("read pw");
    //    [pw pcorr];
        pw_c = [pw combinePWForLoop:ch withCoil:GE_Card_8_new];	// single-ch, real
    //    pw_c = [pw avgForLoop:ch];
        [pw_c saveToFile:@"pw_c" relativePath:YES];
    // for gridding ...
        traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rd, nil];	// traj for 256 x 300
        [traj initRadialTraj];
    //[traj saveAsKOImage:@"pw_test_traj.img"];

        radMap = [traj mapForRadial];   // only dim of traj is used
        [radMap saveAsKOImage:@"pw_test_map.img"];

    TIMER_END("grid");
    } // autorelease pool

    grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
// motion correction loop
	for (i = 0; i < 2; i++) {
        @autoreleasepool {  // inner 2
		printf("==== shift %d =====\n", i);

 	//===== shift corr step ====
    //    [img_c maskWithImage:mask invert:NO smooth:YES];
        [img_c fTriWin3D];  // seems to be OK... no need for manual segmentation

        pwr_c = [img_c reprojectWithMap:radMap];
        [pwr_c copyLoopsOf:pw_c];
        path = [NSString stringWithFormat:@"pwr_%02d.img", i];
        [pwr_c saveToFile:path relativePath:YES];

        switch (mode) {
        case 0 :    // shift
            sft = [pw_c detectShiftWithRef:pwr_c];
            [sft dumpShift:i];
           pws = [pw correctShift:sft];
            break;
        case 1 :    // shift4
            sft = [pw_c detectShift4WithRef:pwr_c];
            [sft dumpShift:i];
            // probably enough to pad/crop within this...
            pws = [pw correctShift4:sft mark:NO over:10];
            break;
        }

        pws_c = [pws combinePWForLoop:ch withCoil:GE_Card_8_new];	// single-ch, real
        [pws_c saveToFile:@"pws" relativePath:YES];

        raw = [pws pwToRaw];
		[raw fft1d:rd direction:REC_INVERSE];   // sin -> raw
        [grid grid2d:raw to:img];
TIMER_END("grid");

        img_c = [img combineForLoop:ch withCoil:GE_Card_8_new];
        path = [NSString stringWithFormat:@"img2_%02d.img", i];
        [img_c saveAsKOImage:path];
        [img_c saveAsRawImage:@"img_c.raw"];	// for Osirix
        } // autorelease pool
    }   // correction loop
    }   // autorelease pool (Outer)

	return 0;
}

