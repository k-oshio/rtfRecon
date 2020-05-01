//
//	rtf3 recon (1 & 2 combined)
//
//	K. Oshio	1-27-2015
//
//  FOV / slthk / nslice is assumed to be 300, 2, 64
//  original pixel size is: xy:1.172, z:2, z/x ratio:1.707
//

/* 2-25-2015

P16384.7	sakai
P10240.7	akasaka

LSB first
RDBM header version = 11.00
 time: 12.735471 (sec) img
====== reproject / correction 0 ======
 time: 19.664232 (sec) pwr
0 1.000000
 time: 11.822875 (sec) img_s
====== reproject / correction 1 ======
 time: 0.723146 (sec) pwr
1 0.219740
 time: 11.768156 (sec) img_s
====== reproject / correction 2 ======
 time: 0.724449 (sec) pwr
2 0.078888
====== scaling ======
 time: 52.026984 (sec) scale/shift
====== final gridding ======
 time: 100.056864 (sec) grid
 time: 218.458848 (sec) TOTAL

*/

#import <RecKit/RecKit.h>
#import <sys/time.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

int
main(int ac, char *av[])
{
	NSString		*path;
	int				pNumber;
    BOOL            zFlip = NO;
	RECCTL			ctl;
	RecImage		*raw, *img, *pw, *mip;
    RecImage        *sft, *sft_prev, *corr;
	RecImage		*dfm;
    RecImage        *img_c, *pw_c, *pwr_c, *pws, *pws_c;
    RecImage        *img_hi, *img_coro;
	RecLoop			*ch, *sl, *pe, *rd;
	RecLoop			*rdZF;
	RecLoop			*kx, *ky;
	RecImage		*radMap, *traj;
    RecGridder      *grid;
    int             i;
    float           dlt, prev_dlt;
    BOOL            term = NO;
    BOOL            scaleBeforeGrid = YES;  // "scale before grid" is better
	int				coil = GE_Card_8_new;
	float			ncc;	// normalized cross correlation

TIMER_ST
    @autoreleasepool {
        system("rm *.img");
        system("rm sft*.txt");
    // raw : ch, sl, pe, rd
        if (ac < 2) {
            printf("rtf3Recon <pNumber(int)><zflip>\n");
            exit(0);
        }
	// (1) === load raw data
        pNumber = atoi(av[1]);
        if (ac > 2) {
            zFlip = YES;
            printf("zFlip\n");
        }
        path = [NSString stringWithFormat:@"P%05d.7", pNumber];
        printf("%s\n", [path UTF8String]);
        raw = [RecImage imageWithPfile:path RECCTL:&ctl];	// 3D
	// raw data loops
        rd = [RecLoop findLoop:@"Read"];
        pe = [RecLoop findLoop:@"Phase"];
        ch = [RecLoop findLoop:@"Channel"];
        sl = [RecLoop findLoop:@"Slice"];
		[raw removePointLoops];

    // rec dim
        rdZF = [RecLoop loopWithName:@"rd_zf" dataLength:256];
        kx   = [RecLoop loopWithName:@"kx"    dataLength:ctl.rc_xres];
        ky   = [RecLoop loopWithName:@"ky"    dataLength:ctl.rc_yres];

        [raw xFlip];    // seq was wrong...
        if (!zFlip) {
            [raw flipForLoop:sl];
        }
	// (2) === pre-processing
    // z-unzip -> handle half fourier (just shift to center, 10240 etc)
        [raw fft1d:sl direction:REC_INVERSE];
    //    [raw shift1d:sl];
		[raw slcSft];
//[raw logP1];
//[raw saveAsKOImage:@"raw.img"];
//exit(0);

    //    [raw fft1d:sl direction:REC_FORWARD];
    // read-axis zero-fill (224 -> 256)
        [raw replaceLoop:rd withLoop:rdZF];
        [raw radPhaseCorr]; // radial phase correction
		[raw saveAsKOImage:@"raw.img"];

//      pw_c = [raw combineForLoop:ch];     // sinogram
//      [pw_c saveAsKOImage:@"pw_sin.img"];

        pw = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, pe, sl, rdZF, nil];
		[raw fft1d:sl direction:REC_FORWARD];
		[raw fft1d:rdZF direction:REC_FORWARD];	// sinogram
        [pw copyImage:raw]; // sinogram -> pw
		pw_c = [pw combinePWForLoop:ch withCoil:coil];
   //     [pw_c saveAsKOImage:@"pw.img"];
        [pw saveAsKOImage:@"pw.img"];

        traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rdZF, nil];
        [traj initRadialTraj];          // for gridding
        radMap = [traj mapForRadial];   // for reprojection

		[traj remove2RR:pw_c thres:0.25];	// reject 2R-R
        [raw fft1d:rdZF direction:REC_INVERSE]; // pw -> raw

	// (3) === initial gridding
        img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
        grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
        [grid grid2d:raw to:img];
        img_c = [img combineForLoop:ch withCoil:coil];
TIMER_END("img");
        [img_c scaleToVal:4000.0];
		[img_c SCIC];
		if (1) {
			[img_c saveAsKOImage:@"img.img"]; // first image (before correction)
			img_hi = img_c;		// not high-res, but used for coro img
		} else {	// hi-res
			img_hi = [img_c scale1dLoop:[img_c xLoop] by:2.0 crop:NO];
			img_hi = [img_hi scale1dLoop:[img_hi yLoop] by:2.0 crop:NO];
			img_hi = [img_hi scale1dLoop:[img_hi zLoop] by:3.4 crop:NO];
			[img_hi saveAsKOImage:@"img.img"]; // first image (before correction)
		}
		img_coro = [RecImage imageWithImage:img_hi];
		[img_coro swapLoop:[img_hi zLoop] withLoop:[img_hi yLoop]];
		[img_coro copyImage:img_hi];
		[img_coro saveAsKOImage:@"img_coro.img"];

    // (4) === iterative shift correction ===
        sft = nil;
        prev_dlt = 2.0;
        for (i = 0; i < 20; i++) {  // max iteration
    //    for (i = 0; i < 1; i++) {
        //    printf("====== reproject / correction %d ======\n", i);
        //    [img_c pwWin];
            [img_c imgWin];
       //     [img_c saveAsKOImage:@"img_w.img"];
            pwr_c = [img_c reprojectWithMap:radMap];	// contains 2dft: img becomes k
            [pwr_c copyLoopsOf:pw_c];
            path = [NSString stringWithFormat:@"pwr%d.img", i];
            [pwr_c saveAsKOImage:path];
            sft_prev = sft;

if (0) {	// 1d ver / hor correlation
	corr = [pw_c xCorrelationWith:pwr_c width:0.2];
	[corr saveAsKOImage:@"img_shearx.img"];
	corr = [pw_c yCorrelationWith:pwr_c width:0.2];
	[corr saveAsKOImage:@"img_sheary.img"];
	exit(0);
}
            corr = [pw_c xyCorrelationWith:pwr_c width:0.2];    // needs tuning
        //    path = [NSString stringWithFormat:@"corr%d.img", i];
            path = @"corr.img";
            [corr saveAsKOImage:path];
       //     sft = [pw_c detectShiftWithCorr:corr];
			sft = [corr corrToSft];
			[sft saveShift:i xDim:[pw xDim] yDim:[pw yDim]];

			[traj removeDeepBreath:sft thres:0.01];
			[grid updateWeight:traj];
		//	grid = [RecGridder gridderWithTrajectory:traj andRecDim:[traj xDim]];


            dlt = [sft deltaWithPreviousSft:sft_prev];
			ncc = [pw_c correlationWith:pwr_c];
            printf("%d %f %f\n", i, dlt, 1.0 - ncc);
			// stopping condition
            if (dlt < 0.1 || dlt > prev_dlt) { // 0.1
                if (scaleBeforeGrid) {
                    break;
                } else {
                    term = YES;
                }
            }
            prev_dlt = dlt;
            pws = [pw correctShift:sft];
            pws_c = [pws combinePWForLoop:ch withCoil:coil];
		//##
		[pws_c gauss2DHP:0.2];

            path = [NSString stringWithFormat:@"pws%d.img", i];
     //       path = @"pws.img";
            [pws_c saveAsKOImage:path];

            raw = [pws pwToSin];
            [raw fft1d:rdZF direction:REC_INVERSE];   // sin -> raw
            [grid grid2d:raw to:img];
[[[grid traj] weight] saveAsKOImage:@"wt.img"];
            img_c = [img combineForLoop:ch withCoil:coil];
            path = [NSString stringWithFormat:@"img_s%d.img", i];
         //   [img_c saveAsKOImage:path];
            [img saveAsKOImage:path];

[img fft2d:REC_INVERSE];
[img logP1];
[img saveAsKOImage:@"img_k.img"];

            if (term) break;
        }
	// (5) === higher order correction
		if (0) {
			printf("dfm correction\n");
			dfm = [pw_c detectDeformWithRef:pwr_c];
			[dfm saveAsKOImage:@"img_dfm.img"];
			pws = [pw correctDeformWithMap:dfm];
			pws_c = [pws combinePWForLoop:ch withCoil:coil];
			[pws_c gauss2DHP:0.2];
			[pws_c saveAsKOImage:@"pws_dfm.img"];
			exit(0);
		}
		// === not done yet ###
	//ncc = [pw_c correlationWith:pwr_c]; printf("ncc = %f\n", -log(1.0 - ncc));

exit(0);
	// (6) === final gridding
        if (scaleBeforeGrid) {
            // scaling
            printf("====== scaling ======\n");
            pw = [pw scale1dLoop:[pw xLoop] by:2.0 crop:NO];
            pw = [pw scale1dLoop:[pw yLoop] by:3.4 crop:NO];
            pws = [pw correctShift:sft];	// ## change to correctDeform later
		//	pws = [pw correctDeform:dfm];
            pws_c = [pws combinePWForLoop:ch withCoil:coil];	// single-ch, real
            [pws_c saveAsKOImage:@"pws_hi.img"];
            corr = [corr correctShift:sft];
            [corr saveAsKOImage:@"corr_s.img"];
TIMER_END("scale/shift");
            raw = [pws pwToSin];
            [raw fft1d:[raw xLoop] direction:REC_INVERSE];   // sin -> raw
            printf("====== final gridding ======\n");
            kx = [RecLoop loopWithDataLength:512];
            ky = [RecLoop loopWithDataLength:512];
            img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, [raw zLoop], ky, kx, nil];
            traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:[raw yLoop], [raw xLoop], nil];
            [traj initRadialTraj];          // for gridding
            grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
            [grid grid2d:raw to:img];
            img_c = [img combineForLoop:ch withCoil:coil];
            [img_c SCIC];
            [img_c saveAsKOImage:@"img_s.img"];
			mip = [img_c partialMipForLoop:[img_c zLoop] depth:16];
			[mip saveAsKOImage:@"img_s_mip.img"];
TIMER_END("grid");
            img_coro = [RecImage imageWithImage:img_c];
            [img_coro swapLoop:[img_c zLoop] withLoop:[img_c yLoop]];
            [img_coro copyImage:img_c];
            [img_coro saveAsKOImage:@"img_s_coro.img"];

            [img_c fft2d:REC_INVERSE];
            [img_c logP1];
            [img_c saveAsKOImage:@"img_k.img"];
        } else {
            // scaling after final gridding
            printf("====== scaling ======\n");
            img_c = [img_c scale1dLoop:[img_c xLoop] by:2.0 crop:NO];
            img_c = [img_c scale1dLoop:[img_c yLoop] by:2.0 crop:NO];
            img_c = [img_c scale1dLoop:[img_c zLoop] by:3.4 crop:NO];
            [img_c SCIC];
            [img_c saveAsKOImage:@"img_s.img"];
            img_coro = [RecImage imageWithImage:img_c];
            [img_coro swapLoop:[img_c zLoop] withLoop:[img_c yLoop]];
            [img_coro copyImage:img_c];
            [img_coro saveAsKOImage:@"img_s_coro.img"];
            pws_c = [pws_c scale1dLoop:[pws_c xLoop] by:2.0 crop:NO];
            pws_c = [pws_c scale1dLoop:[pws_c yLoop] by:3.4 crop:NO];
            [pws_c saveAsKOImage:@"pws_hi.img"];
        }
    } // autoreleasepool
TIMER_TOTAL

	return 0;
}
