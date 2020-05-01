//
//	rtf3 recon 1 (recon without motion correction)
//      rtf3-1 - rtf3-4
//
//	K. Oshio	7-8-2010
//
//  === update
//      step one is just libtest / test2() -> copy
//      add motion correction and update lib
//

#import <RecKit/RecKit.h>
#import <sys/time.h>
#import "../RecKit/timer_macros.h"

int
main(int ac, char *av[])
{
	NSString		*path;
	int				pNumber;
    BOOL            zFlip = NO;
	RecImage		*raw;
	RECCTL			ctl;
	RecImage		*img;
	RecImage		*pw;
	RecLoop			*ch, *sl, *pe, *rd;
	RecLoop			*rdZF;
	RecLoop			*kx, *ky;
	RecImage		*traj;
    RecGridder      *grid;

TIMER_ST
    @autoreleasepool {
    // raw : ch, sl, pe, rd
        if (ac < 2) {
            printf("rtf3Recon <pNumber(int)><zflip>\n");
            exit(0);
        }
        pNumber = atoi(av[1]);
        if (ac > 2) {
            zFlip = YES;
            printf("zFlip\n");
        }
        path = [NSString stringWithFormat:@"P%05d.7", pNumber];
        printf("%s\n", [path UTF8String]);
        raw = [RecImage imageWithPfile:path RECCTL:&ctl];	// 3D
        [raw dumpInfo];
TIMER_END("read raw");
        [raw saveToFile:@"raw" relativePath:YES];
TIMER_END("save raw");

        rd = [RecLoop findLoop:@"Read"];
        pe = [RecLoop findLoop:@"Phase"];
        ch = [RecLoop findLoop:@"Channel"];
        sl = [RecLoop findLoop:@"Slice"];

    // rec dim
        rdZF = [RecLoop loopWithName:@"rd_zf" dataLength:256];
        kx = [RecLoop loopWithName:@"kx" dataLength:ctl.rc_xres];
        ky = [RecLoop loopWithName:@"ky" dataLength:ctl.rc_yres];

        [raw xFlip];
        // seq was wrong...
        if (!zFlip) {
            [raw flipForLoop:sl];
        }

        [raw fft1d:sl direction:REC_INVERSE];   // img -> raw
        [raw shift1d:sl];
TIMER_END("slice IFT");
        [raw fft1d:sl direction:REC_FORWARD];   // raw -> img
TIMER_END("slice FT");
        [raw replaceLoop:rd withLoop:rdZF];
        [raw radPhaseCorr]; // includes read-FT
TIMER_END("rad pcorr");
//        pw = [raw combineForLoop:ch];     // sinogram
//        [pw saveToFile:@"pw" relativePath:YES];

        pw = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, pe, sl, rdZF, nil];
        [pw copyImage:raw];
        [pw saveToFile:@"pw" relativePath:YES];
        pw = [pw combinePWForLoop:ch withCoil:GE_Card_8_new];
        [pw saveToFile:@"pwc" relativePath:YES];
TIMER_END("save PW");
        traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rdZF, nil];	// traj for 256 x 300
        [traj initRadialTraj];
        [raw fft1d:rdZF direction:REC_INVERSE]; // img -> raw

        img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
        grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
        [grid grid2d:raw to:img];
//        [grid grid2d:raw to:img iteration:2];
TIMER_END("grid");
        [img saveToFile:@"img" relativePath:YES];	// RecImage obj archive7

        img = [img combineForLoop:ch withCoil:GE_Card_8_new];
    //    [img pcorr];
    //    img = [img complexCombineForLoop:ch withCoil:GE_Card_8_new];
TIMER_END("save img");
        [img scaleToVal:4000.0];
        [img saveAsKOImage:@"img.img"];
        [img saveAsRawImage:@"img.raw"];	// for Osirix
    } /// autoreleasepool
TIMER_TOTAL

	return 0;
}
