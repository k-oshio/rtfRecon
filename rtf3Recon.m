//
//	rtf3 recon (main, GE)
//
//	K. Oshio	1-27-2015
//
//  FOV / slthk / nslice is assumed to be 300, 2, 64
//  original pixel size is: xy:1.172, z:2, z/x ratio:1.707
//
//  forked off from rtf3Recon4.m    4-28-2020

//  started update for central line nav 5-2-2020
//  cmd-line arg -> hard-coded path + zflip etc


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
	RecImage		*raw, *img, *pw, *mask;
    RecImage        *sft, *sft_prev, *corr;
	RecImage		*dfm;
    RecImage        *img_c, *pw_c, *pwr, *pws, *pws_c;
    RecImage        *img_hi, *img_coro, *img_coro_c, *img_lap, *map_coro;
	RecLoop			*ch, *sl, *pe, *rd;
	RecLoop			*rdZF;
	RecLoop			*kx, *ky;
	RecImage		*radMap, *traj, *thTab;
    RecImage        *imgs;
	RecImage		*tmp_img;
    RecGridder      *grid;
    int             i;
    float           dlt, prev_dlt;
	float			sft_scl;
	int				coil = GE_Card_8_new;
 
TIMER_ST
    @autoreleasepool {
        if (0) {
            pw = [RecImage imageFromFile:@"pw_sav.recimg" relativePath:YES];
            [pw dumpLoops];
            sft = [pw shiftFromK0];
            exit(0);
        }
        
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
        [raw shift1d:sl];

    // read-axis zero-fill (224 -> 256)
        [raw replaceLoop:rd withLoop:rdZF];
        [raw radPhaseCorr]; // radial phase correction
		[raw saveAsKOImage:@"raw.img"];

        pw = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, pe, sl, rdZF, nil];

		[raw fft1d:sl direction:REC_FORWARD];
		[raw fft1d:rdZF direction:REC_FORWARD];	// sinogram
        [pw copyImage:raw]; // sinogram -> pw
		pw_c = [pw combinePWForLoop:ch withCoil:coil];
        [pw_c saveAsKOImage:@"pw_c.img"];
        [pw saveAsKOImage:@"pw.img"];
        [pw saveToFile:@"pw_sav" relativePath:YES];
   
        traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rdZF, nil];

        thTab = [traj initRadialTraj];          // for gridding
        radMap = [traj mapForRadial];   // for reprojection
        [raw fft1d:rdZF direction:REC_INVERSE]; // pw -> raw

	// (3) === initial gridding
        img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
        grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
        [grid grid2d:raw to:img];
        img_c = [img combineForLoop:ch withCoil:coil];
   //     [img_c scaleToVal:4000.0];
	//	[img_c SCIC];
		[img_c saveAsKOImage:@"img_c.img"]; // first image (before correction)

	//	mask = [img_c copy];
	//	[mask gauss3DLP:0.1];
	//	[mask saveAsKOImage:@"mask_1.img"];

		img_coro_c = [RecImage imageWithImage:img_c];
		[img_coro_c swapLoop:[img_c zLoop] withLoop:[img_c yLoop]];
		[img_coro_c copyImage:img_c];
		[img_coro_c saveAsKOImage:@"img_coro.img"];

    // sft est
        sft = [pw shiftFromK0];
        [sft saveAsKOImage:@"IMG_sft"];
    
    // correction (just copied)
    // -> move to lib (category)
    // input: sft, pw, grid     output:imgs
   //     - (RecImage *)stepCorrWithPW:(RecImage *)pw gridder:(RecGridder *)grid sft:(RecImage *)sft nScale:(int)nScl // input(self) is img
        imgs = [img stepCorrWithPW:pw gridder:grid sft:sft];
        [imgs saveAsKOImage:@"IMG_imgs"];
/*
        if (1) {
            float            scl;
            int                nscl = 10;
            float            scl_range = 3.0;    // (0 - 2.0)
            RecImage        *imgs, *mxSft;
            RecImage        *sft_scl, *pws;
            RecLoop            *scLp;
            float            lap_w = 0.2;    // 0.3

            printf("4: step-shift\n");
            grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim] nop:1 densityCorrection:YES];

            scLp = [RecLoop loopWithDataLength:nscl];
            imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];    // coro
        //    imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:[img yLoop], scLp, [img zLoop], [img xLoop], nil];    // coro, scale first (doesn't work)

            for (i = 0; i < nscl; i++) {
                scl = (float)i * scl_range / nscl;
                printf("scale %d = %f\n", i, scl);

                sft_scl = [sft copy];            // pixels
                [sft_scl multByConst:scl];
                pws = [pw correctZShift:sft_scl]; // z shift

                tmp_img = [pws copy];
                raw = [pws copy];
                [raw fft1d:[raw xLoop] direction:REC_INVERSE];
                [raw swapLoop:[raw yLoop] withLoop:[raw zLoop]];
                [grid grid2d:raw to:img];
                img_c = [img combineForLoop:ch];
                [imgs copySlice:img_c atIndex:i forLoop:scLp];
            }
            tmp_img = [imgs copy];
            [tmp_img swapLoop:scLp withLoop:[img yLoop]];
            [tmp_img saveAsKOImage:@"imgs_scl.sav"];

// (5) === find best focus scale ====
            printf("5: select best shift\n");
            tmp_img = [imgs copy];
//            [tmp_img laplace2d:REC_FORWARD];
            [tmp_img grad1dForLoop:[tmp_img yLoop]];
            [tmp_img magnitude];
            [tmp_img square];
            [tmp_img gauss3DLP:lap_w];
            [tmp_img swapLoop:scLp withLoop:[img yLoop]];    // ## chk
            [tmp_img saveAsKOImage:@"imgs_scl_lap.img"];

            mxSft = [tmp_img peakIndexForLoop:scLp];    // sharper with larger laplacian
            [mxSft saveAsKOImage:@"imgMxSft.img"];        // ok
            [imgs swapLoop:scLp withLoop:[img yLoop]];    // ## chk
            tmp_img = [imgs selectSft:mxSft];            // ### loop ? [scl z x]
            // ####
            [tmp_img saveAsKOImage:@"imgs_final.img"];

            [mxSft gauss2DLP:0.1];
            [mxSft saveAsKOImage:@"imgMxSft_f.img"];
            tmp_img = [imgs selectSftF:mxSft];
            [tmp_img saveAsKOImage:@"imgs_final_F.img"];
            [tmp_img swapLoop:[tmp_img yLoop] withLoop:[tmp_img zLoop]];
            [tmp_img saveAsKOImage:@"imgs_final_ax.img"];
        } // === focused rigid-body correction === 
*/


    } // autoreleasepool
TIMER_TOTAL

	return 0;
}
