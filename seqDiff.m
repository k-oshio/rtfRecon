//
//  quick test routines
//

#import <RecKit/RecKit.h>
#import "RecImagePW.h"

int
main()
{
    @autoreleasepool {

    // slice nopw using coil profile
        if (0) {
            RecImage    *pw = [RecImage imageFromFile:@"pw" relativePath:YES];
            RecImage    *pw_e = [pw expandSliceBy:10];
            [pw_e saveAsKOImage:@"pw_e.img"];
            [pw_e cropBy:10];
            [pw_e saveAsKOImage:@"pw_ee.img"];
        }
    // filter before combining -> doesn't work
        if (0) {
            RecImage    *pw = [RecImage imageFromFile:@"pw" relativePath:YES];
            RecImage    *pwc_f;
            RecLoop     *ch = [RecLoop findLoop:@"Channel"];

            [pw magnitude];
            [pw rhoFilter];
            [pw makeComplex];
            pwc_f = [pw combinePWForLoop:ch withCoil:GE_Card_8_new];
            [pwc_f saveAsKOImage:@"pwc_f.img"];
        }

    // HP-filter along proj (non-pow-of-2)
        if (0) {
            RecImage    *sin_c, *sin_v;
            RecLoop     *ky, *ky_2;

            sin_c = [RecImage imageFromFile:@"sin_c" relativePath:YES];
            // zero-fill
            ky_2 = [RecLoop loopWithDataLength:512];
            ky = [sin_c yLoop];
            [sin_c replaceLoop:ky withLoop:ky_2];
            // HPF
            [sin_c gauss2DHP:1.0];
            // crop
            [sin_c replaceLoop:ky_2 withLoop:ky];
            [sin_c saveAsKOImage:@"sin_f.img"];
            // measure variance along proj
        //    sin_v = [sin_c detectMotion];
            [sin_v saveAsKOImage:@"sin_v.img"];
        }

        return 0;
    }
}

