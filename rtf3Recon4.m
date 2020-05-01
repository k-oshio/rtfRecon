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
	RecImage		*raw, *img, *pw, *mask;
    RecImage        *sft, *sft_prev, *corr;
	RecImage		*dfm;
    RecImage        *img_c, *pw_c, *pwr, *pws, *pws_c;
    RecImage        *img_hi, *img_coro, *img_coro_c, *img_lap, *map_coro;
	RecLoop			*ch, *sl, *pe, *rd;
	RecLoop			*rdZF;
	RecLoop			*kx, *ky;
	RecImage		*radMap, *traj, *thTab;
	RecImage		*tmp_img;
    RecGridder      *grid;
    int             i;
    float           dlt, prev_dlt;
	float			sft_scl;
    BOOL            term = NO;
    BOOL            scaleBeforeGrid = YES;  // "scale before grid" is better
	BOOL			nav = NO;
	int				coil = GE_Card_8_new;
	float			ncc;	// normalized cross correlation
	RecBlock3		blk;

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


exit(0);
// 4-28-2020


	// virtual navigator ###### -> output is [sft]
	// scale: diaphragm motion is largest
		if (nav) {	// left diaphragm (spleen) -> optical flow
			RecImage	*nav, *ref, *srt;
			int			i, n;
			float		*p, mx;
			
	ref = [pw sliceAtIndex:4 forLoop:ch];
	[ref saveAsKOImage:@"ref.img"];
			nav = [ref avgForLoop:[ref xLoop]];
			[nav saveAsKOImage:@"nav.img"];
			
	//		sft = [nav opticalFlow1d];
			sft = [nav detectShift];
			[sft saveAsKOImage:@"sft.img"];
			if (0) {
				p = [sft data];
				n = [sft dataLength];
				for (i = 0; i < n; i++) {
					printf("%d %f\n", i, p[i]);
				}
			}

			ref = [ref correctZShift:sft];	// z
			[ref saveAsKOImage:@"pw_s.img"];

			srt = [sft sortZSft];
			[srt saveAsKOImage:@"sft_sort.img"];
			if (0) {
				p = [srt data];
				n = [srt dataLength];
				for (i = 0; i < n; i++) {
					printf("%d %f\n", i, p[i]);
				}
			}
			sft_scl = 1.2;
		} else {
	// PCA shift detection ######
	// scale estimation is necessary
	// chk direction
			pwr = [img reprojectWithMap:radMap];
			corr = [pw xyCorrelationWith:pwr width:0.2 triFilt:YES];	// 0.2
			sft = [corr corrToSftZ];	// pixels
		//	[sft negate];	// ??
			[sft saveAsKOImage:@"sft.img"];
			[sft dumpShift];
		// ### testing...
			sft_scl = 1.0;	// ### est scale & direction using binned recon
			[sft multByConst:sft_scl];
			pws = [pw_c correctZShift:sft];
			[pws saveAsKOImage:@"pws.img"];
		}

		img_lap = [img_coro_c copy];
		[img_lap laplace3d:REC_FORWARD];
		[img_lap magnitude];

		[img_lap gauss3DLP:0.1];
		[img_lap saveAsKOImage:@"img_coro_l.img"];

	// binned recon / 1D non-rigid
		if (1) {
			int			nBin = 2;
			RecImage	*pwBins, *bTab, *map;
			RecImage	*pwb, *rawb, *trajb;	// probably no need for pwb
			RecImage	*imgb[3], *def1, *def2, *img1, *img2, *img0;
			RecImage	*sftb;
			int			i, j;
			float		mn;

			// pwBins [proj]
			pwBins = [sft binTabWithNBins:nBin binSize:[sft xDim] / 2];
			[pwBins saveAsKOImage:@"bins.img"];

			// reconstruct each bin
			for (i = 0; i < nBin; i++) {
				bTab = [pwBins sliceAtIndex:i forLoop:[pwBins yLoop]];
				pwb = [pw replaceLoop:pe withTab:bTab];	// 0: exp, 1: mid, 2 inh
				rawb = [pwb copy];
				[rawb fft1d:[rawb xLoop] direction:REC_INVERSE];
				[rawb swapLoop:sl withLoop:[pwb zLoop]];

				trajb = [traj replaceLoop:pe withTab:bTab];
				grid = [RecGridder gridderWithTrajectory:trajb andRecDim:[img xDim]];
				[grid grid2d:rawb to:img];
				path = [NSString stringWithFormat:@"imgb%d.img", i];
			//	[img saveAsKOImage:path];
				[img saveComb:ch as:path];
				imgb[i] = [img copy];	// save

				img_coro_c = [img combineForLoop:ch];
				[img_coro_c swapLoop:[img_coro_c zLoop] withLoop:[img_coro_c yLoop]];
				path = [NSString stringWithFormat:@"imgb%d_coro.img", i];
				[img_coro_c saveAsKOImage:path];
			//##
				[img_coro_c laplace3d:REC_FORWARD];
				[img_coro_c magnitude];
	//	[img_coro_c gradMag3d];
				[img_coro_c gauss3DLP:0.1];
				[img_coro_c divImage:img_lap];
			//	[img_coro_c subImage:img_lap];
				path = [NSString stringWithFormat:@"imgb%d_coro_l.img", i];
				[img_coro_c saveAsKOImage:path];
			}

	// 2 step diff
			if (1) {
				float		*df, *dfz;
				float		mxz;
				int			i, n;
				RecImage	*sftz, *sftxy, *hist;

				// correction (just add 3 bins after warping)
				// warping channels separately doesn't help
				img1 = [imgb[1] copy];
				img0 = [imgb[0] copy];	// ref
				img0 = [img0 combineForLoop:ch];
				img1 = [img1 combineForLoop:ch];
				def1 = [img1 defVectorWithRef:img0];
				[def1 saveAsKOImage:@"def1.img"];
				[def1 negate];

				map = [def1 scaleMapToImage:img1];
				map_coro = [map copy];
				[map_coro swapLoop:[map_coro yLoop] withLoop:[map_coro zLoop]];
				[map_coro saveAsKOImage:@"def10.img"];

				// ### make 2d histo (zy, zx)
				sftz = [RecImage imageOfType:RECIMAGE_REAL withImage:map_coro];
				n = [sftz dataLength];
				df = [sftz data];
				dfz = [map_coro data] + n * 2;	// z
				for (i = 0; i < n; i++) {
					df[i] = dfz[i];
				}
				sftxy = [RecImage imageWithImage:sftz];
				df = [sftxy data];
				dfz = [map_coro data] + n;
				for (i = 0; i < n; i++) {
					df[i] = dfz[i];
				}
				[sftz saveAsKOImage:@"sftz.img"];
				[sftxy saveAsKOImage:@"sftx.img"];
				hist = [RecImage imageOfType:RECIMAGE_REAL
					withLoops:[sftz zLoop], [RecLoop loopWithDim:100], [RecLoop loopWithDim:100], nil];
				[hist histogram2dWithX:sftz andY:sftxy];
			//	hist = [hist logP1];
				[hist saveAsKOImage:@"sft_hist.img"];
exit(0);

				map = dispToMap(map);
				tmp_img = [RecImage imageWithImage:img0];
				[tmp_img resample3d:img1 withMap:map];
				
				[tmp_img saveAsKOImage:@"imgs10.img"];
				[img_coro_c copyImage:tmp_img];
				[img_coro_c saveAsKOImage:@"imgs10_coro.img"];
				
				n = [def1 dataLength];
				df = [def1 data] + n * 2;	// z shift
				mxz = 0;
				for (i = 0; i < n; i++) {
					if (fabs(mxz) < fabs(df[i])) {
						mxz = df[i];
					}
				}
				printf("mx zshift = %f\n", mxz);
			}

		// focused rigid body corr
			if (1) {
				float			scl;
				int				nscl = 10;	// 20
				RecImage		*fcs, *imgs, *mxSft;
				RecLoop			*scLp;
				RecLoopControl	*lc;

				scLp = [RecLoop loopWithDataLength:nscl];
				// coronal
				fcs  = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];
				imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];

				grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
				for (i = 0; i < nscl; i++) {
					scl = sft_scl / nscl * (i + 1);
					printf("scale %d = %f\n", i, scl);
					sftb = [sft copy];
					[sftb multByConst:scl];
					pws = [pw correctZShift:sftb];
					tmp_img = [pws copy];
					tmp_img = [tmp_img combineForLoop:ch];
					path = [NSString stringWithFormat:@"pws%d.img", i];
					[tmp_img saveAsKOImage:path];

					raw = [pws copy];
					[raw fft1d:[raw xLoop] direction:REC_INVERSE];
					[raw swapLoop:[raw yLoop] withLoop:[raw zLoop]];

					img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
					[grid grid2d:raw to:img];
					img_c = [img combineForLoop:ch withCoil:coil];
					path = [NSString stringWithFormat:@"imgs%d.img", i];
					[img_c saveAsKOImage:path];
					[imgs copySlice:img_c atIndex:i forLoop:scLp];
					img_coro_c = [img_c copy];

					[img_coro_c laplace3d:REC_FORWARD];
					[img_coro_c magnitude];
					[img_coro_c gauss3DLP:0.1];
					[img_coro_c divImage:img_lap];
					[fcs copySlice:img_coro_c atIndex:i forLoop:scLp];
				}
				[fcs saveAsKOImage:@"imgLapScl.img"];
				[imgs saveAsKOImage:@"imgs.img"];

			// find best focus index
				mxSft = [fcs maxIndexForLoop:scLp];
				[mxSft saveAsKOImage:@"imgMxSft.img"];

			// select best image
				img_coro = [imgs selectSft:mxSft];
				[img_coro saveAsKOImage:@"img_final_coro.img"];
				[img_c copyImage:img_coro];
				[img_c saveAsKOImage:@"img_final.img"];
			}
			exit(0);
		}







    } // autoreleasepool
TIMER_TOTAL

	return 0;
}
