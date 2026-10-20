#!/usr/bin/env python3
"""
reproduce.py - companion to the research note "The recycled numbers inside Certicom's curve generators".

Python 3.8+, standard library only. Runs in about a minute.

What it does
  For every standard curve below, it computes every "hidden point" P with k*P = G, where G is the published generator
  and k is a divisor of the cofactor h greater than 1 (k = 2 for curves of prime order, i.e. "half the generator").
  On curves with h > 1 this includes points outside the prime-order subgroup. Then it runs four tests:
    1. Which hidden x-coordinates are far shorter than chance allows?            (known before this note; credited)
    2. Blind, all pairs of curves: which hidden x-coordinates share a long run of hex digits?
    3. Blind: which hidden x-coordinates share a long run with a PUBLISHED generator x of another curve?
    4. Blind: which hidden x-coordinates (and published generators) have hex digits far too uneven for chance?
  and finally checks each specific claim made in the note.

Curve parameters: copied from OpenSSL's crypto/ec/ec_curve.c (master, October 2026), which reproduces SEC 2 v1.0 (2000),
ANSI X9.62 (1998) and FIPS 186. Every generator is checked to lie on its curve and to have order n before use.
The six X9.62 curves with cofactors in the thousands (c2pnb176w1, c2pnb208w1, c2pnb272w1, c2pnb304w1, c2pnb368w1,
c2tnb431r1) are left out of the blind tests for speed (the note reports the full search of them separately); one specific
c2tnb431r1 point is checked directly at the end. c2tnb359v1 (cofactor 76) is included.
"""
from itertools import combinations
import math

# (name, field ('p' prime / '2' binary), cofactor h, p or reduction polynomial, a, b, Gx, Gy, n)
CURVES = [
    ('secp192r1 (P-192)', 'p', 1, 0xfffffffffffffffffffffffffffffffeffffffffffffffff, 0xfffffffffffffffffffffffffffffffefffffffffffffffc, 0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1,
     0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012,
     0x7192b95ffc8da78631011ed6b24cdd573f977a11e794811,
     0xffffffffffffffffffffffff99def836146bc9b1b4d22831),
    ('secp224r1 (P-224)', 'p', 1, 0xffffffffffffffffffffffffffffffff000000000000000000000001, 0xfffffffffffffffffffffffffffffffefffffffffffffffffffffffe, 0xb4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4,
     0xb70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21,
     0xbd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34,
     0xffffffffffffffffffffffffffff16a2e0b8f03e13dd29455c5c2a3d),
    ('secp384r1 (P-384)', 'p', 1, 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff, 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffc, 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef,
     0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7,
     0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f,
     0xffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973),
    ('secp521r1 (P-521)', 'p', 1, 0x1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff, 0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc, 0x51953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00,
     0xc6858e06b70404e9cd9e3ecb662395b4429c648139053fb521f828af606b4d3dbaa14b5e77efe75928fe1dc127a2ffa8de3348b3c1856a429bf97e7e31c2e5bd66,
     0x11839296a789a3bc0045c8a5fb42c7d1bd998f54449579b446817afbd17273e662c97ee72995ef42640c550b9013fad0761353c7086a272c24088be94769fd16650,
     0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa51868783bf2f966b7fcc0148f709a5d03bb5c9b8899c47aebb6fb71e91386409),
    ('prime192v2', 'p', 1, 0xfffffffffffffffffffffffffffffffeffffffffffffffff, 0xfffffffffffffffffffffffffffffffefffffffffffffffc, 0xcc22d6dfb95c6b25e49c0d6364a4e5980c393aa21668d953,
     0xeea2bae7e1497842f2de7769cfe9c989c072ad696f48034a,
     0x6574d11d69b6ec7a672bb82a083df2f2b0847de970b2de15,
     0xfffffffffffffffffffffffe5fb1a724dc80418648d8dd31),
    ('prime192v3', 'p', 1, 0xfffffffffffffffffffffffffffffffeffffffffffffffff, 0xfffffffffffffffffffffffffffffffefffffffffffffffc, 0x22123dc2395a05caa7423daeccc94760a7d462256bd56916,
     0x7d29778100c65a1da1783716588dce2b8b4aee8e228f1896,
     0x38a90f22637337334b49dcb66a6dc8f9978aca7648a943b0,
     0xffffffffffffffffffffffff7a62d031c83f4294f640ec13),
    ('prime239v1', 'p', 1, 0x7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff, 0x7fffffffffffffffffffffff7fffffffffff8000000000007ffffffffffc, 0x6b016c3bdcf18941d0d654921475ca71a9db2fb27d1d37796185c2942c0a,
     0xffa963cdca8816ccc33b8642bedf905c3d358573d3f27fbbd3b3cb9aaaf,
     0x7debe8e4e90a5dae6e4054ca530ba04654b36818ce226b39fccb7b02f1ae,
     0x7fffffffffffffffffffffff7fffff9e5e9a9f5d9071fbd1522688909d0b),
    ('prime239v2', 'p', 1, 0x7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff, 0x7fffffffffffffffffffffff7fffffffffff8000000000007ffffffffffc, 0x617fab6832576cbbfed50d99f0249c3fee58b94ba0038c7ae84c8c832f2c,
     0x38af09d98727705120c921bb5e9e26296a3cdcf2f35757a0eafd87b830e7,
     0x5b0125e4dbea0ec7206da0fc01d9b081329fb555de6ef460237dff8be4ba,
     0x7fffffffffffffffffffffff800000cfa7e8594377d414c03821bc582063),
    ('prime239v3', 'p', 1, 0x7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff, 0x7fffffffffffffffffffffff7fffffffffff8000000000007ffffffffffc, 0x255705fa2a306654b1f4cb03d6a750a30c250102d4988717d9ba15ab6d3e,
     0x6768ae8e18bb92cfcf005c949aa2c6d94853d0e660bbf854b1c9505fe95a,
     0x1607e6898f390c06bc1d552bad226f3b6fcfe48b6e818499af18e3ed6cf3,
     0x7fffffffffffffffffffffff7fffff975deb41b3a6057c3c432146526551),
    ('secp256r1 (P-256)', 'p', 1, 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff, 0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc, 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
     0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
     0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5,
     0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551),
    ('secp112r1', 'p', 1, 0xdb7c2abf62e35e668076bead208b, 0xdb7c2abf62e35e668076bead2088, 0x659ef8ba043916eede8911702b22,
     0x9487239995a5ee76b55f9c2f098,
     0xa89ce5af8724c0a23e0e0ff77500,
     0xdb7c2abf62e35e7628dfac6561c5),
    ('secp112r2', 'p', 4, 0xdb7c2abf62e35e668076bead208b, 0x6127c24c05f38a0aaaf65c0ef02c, 0x51def1815db5ed74fcc34c85d709,
     0x4ba30ab5e892b4e1649dd0928643,
     0xadcd46f5882e3747def36e956e97,
     0x36df0aafd8b8d7597ca10520d04b),
    ('secp128r1', 'p', 1, 0xfffffffdffffffffffffffffffffffff, 0xfffffffdfffffffffffffffffffffffc, 0xe87579c11079f43dd824993c2cee5ed3,
     0x161ff7528b899b2d0c28607ca52c5b86,
     0xcf5ac8395bafeb13c02da292dded7a83,
     0xfffffffe0000000075a30d1b9038a115),
    ('secp128r2', 'p', 4, 0xfffffffdffffffffffffffffffffffff, 0xd6031998d1b3bbfebf59cc9bbff9aee1, 0x5eeefca380d02919dc2c6558bb6d8a5d,
     0x7b6aa5d85e572983e6fb32a7cdebc140,
     0x27b6916a894d3aee7106fe805fc34b44,
     0x3fffffff7fffffffbe0024720613b5a3),
    ('secp160k1', 'p', 1, 0xfffffffffffffffffffffffffffffffeffffac73, 0x0, 0x7,
     0x3b4c382ce37aa192a4019e763036f4f5dd4d7ebb,
     0x938cf935318fdced6bc28286531733c3f03c4fee,
     0x100000000000000000001b8fa16dfab9aca16b6b3),
    ('secp160r1', 'p', 1, 0xffffffffffffffffffffffffffffffff7fffffff, 0xffffffffffffffffffffffffffffffff7ffffffc, 0x1c97befc54bd7a8b65acf89f81d4d4adc565fa45,
     0x4a96b5688ef573284664698968c38bb913cbfc82,
     0x23a628553168947d59dcc912042351377ac5fb32,
     0x100000000000000000001f4c8f927aed3ca752257),
    ('secp160r2', 'p', 1, 0xfffffffffffffffffffffffffffffffeffffac73, 0xfffffffffffffffffffffffffffffffeffffac70, 0xb4e134d3fb59eb8bab57274904664d5af50388ba,
     0x52dcb034293a117e1f4ff11b30f7199d3144ce6d,
     0xfeaffef2e331f296e071fa0df9982cfea7d43f2e,
     0x100000000000000000000351ee786a818f3a1a16b),
    ('secp192k1', 'p', 1, 0xfffffffffffffffffffffffffffffffffffffffeffffee37, 0x0, 0x3,
     0xdb4ff10ec057e9ae26b07d0280b7f4341da5d1b1eae06c7d,
     0x9b2f2f6d9c5628a7844163d015be86344082aa88d95e2f9d,
     0xfffffffffffffffffffffffe26f2fc170f69466a74defd8d),
    ('secp224k1', 'p', 1, 0xfffffffffffffffffffffffffffffffffffffffffffffffeffffe56d, 0x0, 0x5,
     0xa1455b334df099df30fc28a169a467e9e47075a90f7e650eb6b7a45c,
     0x7e089fed7fba344282cafbd6f7e319f7c0b0bd59e2ca4bdb556d61a5,
     0x10000000000000000000000000001dce8d2ec6184caf0a971769fb1f7),
    ('secp256k1', 'p', 1, 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f, 0x0, 0x7,
     0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
     0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8,
     0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141),
    ('sect113r1', '2', 2, 0x20000000000000000000000000201, 0x3088250ca6e7c7fe649ce85820f7, 0xe8bee4d3e2260744188be0e9c723,
     0x9d73616f35f4ab1407d73562c10f,
     0xa52830277958ee84d1315ed31886,
     0x100000000000000d9ccec8a39e56f),
    ('sect113r2', '2', 2, 0x20000000000000000000000000201, 0x689918dbec7e5a0dd6dfc0aa55c7, 0x95e9a9ec9b297bd4bf36e059184f,
     0x1a57a6a7b26ca5ef52fcdb8164797,
     0xb3adc94ed1fe674c06e695baba1d,
     0x10000000000000108789b2496af93),
    ('sect131r1', '2', 2, 0x80000000000000000000000000000010d, 0x7a11b09a76b562144418ff3ff8c2570b8, 0x217c05610884b63b9c6c7291678f9d341,
     0x81baf91fdf9833c40f9c181343638399,
     0x78c6e7ea38c001f73c8134b1b4ef9e150,
     0x400000000000000023123953a9464b54d),
    ('sect131r2', '2', 2, 0x80000000000000000000000000000010d, 0x3e5a88919d7cafcbf415f07c2176573b2, 0x4b8266a46c55657ac734ce38f018f2192,
     0x356dcd8f2f95031ad652d23951bb366a8,
     0x648f06d867940a5366d9e265de9eb240f,
     0x400000000000000016954a233049ba98f),
    ('sect163k1 (K-163)', '2', 2, 0x800000000000000000000000000000000000000c9, 0x1, 0x1,
     0x2fe13c0537bbc11acaa07d793de4e6d5e5c94eee8,
     0x289070fb05d38ff58321f2e800536d538ccdaa3d9,
     0x4000000000000000000020108a2e0cc0d99f8a5ef),
    ('sect163r1', '2', 2, 0x800000000000000000000000000000000000000c9, 0x7b6882caaefa84f9554ff8428bd88e246d2782ae2, 0x713612dcddcb40aab946bda29ca91f73af958afd9,
     0x369979697ab43897789566789567f787a7876a654,
     0x435edb42efafb2989d51fefce3c80988f41ff883,
     0x3ffffffffffffffffffff48aab689c29ca710279b),
    ('sect163r2 (B-163)', '2', 2, 0x800000000000000000000000000000000000000c9, 0x1, 0x20a601907b8c953ca1481eb10512f78744a3205fd,
     0x3f0eba16286a2d57ea0991168d4994637e8343e36,
     0xd51fbc6c71a0094fa2cdd545b11c5c0c797324f1,
     0x40000000000000000000292fe77e70c12a4234c33),
    ('sect193r1', '2', 2, 0x2000000000000000000000000000000000000000000008001, 0x17858feb7a98975169e171f77b4087de098ac8a911df7b01, 0xfdfb49bfe6c3a89facadaa7a1e5bbc7cc1c2e5d831478814,
     0x1f481bc5f0ff84a74ad6cdf6fdef4bf6179625372d8c0c5e1,
     0x25e399f2903712ccf3ea9e3a1ad17fb0b3201b6af7ce1b05,
     0x1000000000000000000000000c7f34a778f443acc920eba49),
    ('sect193r2', '2', 2, 0x2000000000000000000000000000000000000000000008001, 0x163f35a5137c2ce3ea6ed8667190b0bc43ecd69977702709b, 0xc9bb9e8927d4d64c377e2ab2856a5b16e3efb7f61d4316ae,
     0xd9b67d192e0367c803f39e1a7e82ca14a651350aae617e8f,
     0x1ce94335607c304ac29e7defbd9ca01f596f927224cdecf6c,
     0x10000000000000000000000015aab561b005413ccd4ee99d5),
    ('sect233k1 (K-233)', '2', 4, 0x20000000000000000000000000000000000000004000000000000000001, 0x0, 0x1,
     0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,
     0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3,
     0x8000000000000000000000000000069d5bb915bcd46efb1ad5f173abdf),
    ('sect233r1 (B-233)', '2', 2, 0x20000000000000000000000000000000000000004000000000000000001, 0x1, 0x66647ede6c332c7f8c0923bb58213b333b20e9ce4281fe115f7d8f90ad,
     0xfac9dfcbac8313bb2139f1bb755fef65bc391f8b36f8f8eb7371fd558b,
     0x1006a08a41903350678e58528bebf8a0beff867a7ca36716f7e01f81052,
     0x1000000000000000000000000000013e974e72f8a6922031d2603cfe0d7),
    ('sect239k1', '2', 4, 0x800000000000000000004000000000000000000000000000000000000001, 0x0, 0x1,
     0x29a0b6a887a983e9730988a68727a8b2d126c44cc2cc7b2a6555193035dc,
     0x76310804f12e549bdb011c103089e73510acb275fc312a5dc6b76553f0ca,
     0x2000000000000000000000000000005a79fec67cb6e91f1c1da800e478a5),
    ('sect283k1 (K-283)', '2', 4, 0x800000000000000000000000000000000000000000000000000000000000000000010a1, 0x0, 0x1,
     0x503213f78ca44883f1a3b8162f188e553cd265f23c1567a16876913b0c2ac2458492836,
     0x1ccda380f1c9e318d90f95d07e5426fe87e45c0e8184698e45962364e34116177dd2259,
     0x1ffffffffffffffffffffffffffffffffffe9ae2ed07577265dff7f94451e061e163c61),
    ('sect283r1 (B-283)', '2', 2, 0x800000000000000000000000000000000000000000000000000000000000000000010a1, 0x1, 0x27b680ac8b8596da5a4af8a19a0303fca97fd7645309fa2a581485af6263e313b79a2f5,
     0x5f939258db7dd90e1934f8c70b0dfec2eed25b8557eac9c80e2e198f8cdbecd86b12053,
     0x3676854fe24141cb98fe6d4b20d02b4516ff702350eddb0826779c813f0df45be8112f4,
     0x3ffffffffffffffffffffffffffffffffffef90399660fc938a90165b042a7cefadb307),
    ('sect409k1 (K-409)', '2', 4, 0x2000000000000000000000000000000000000000000000000000000000000000000000000000000008000000000000000000001, 0x0, 0x1,
     0x60f05f658f49c1ad3ab1890f7184210efd0987e307c84c27accfb8f9f67cc2c460189eb5aaaa62ee222eb1b35540cfe9023746,
     0x1e369050b7c4e42acba1dacbf04299c3460782f918ea427e6325165e9ea10e3da5f6c42e9c55215aa9ca27a5863ec48d8e0286b,
     0x7ffffffffffffffffffffffffffffffffffffffffffffffffffe5f83b2d4ea20400ec4557d5ed3e3e7ca5b4b5c83b8e01e5fcf),
    ('sect409r1 (B-409)', '2', 2, 0x2000000000000000000000000000000000000000000000000000000000000000000000000000000008000000000000000000001, 0x1, 0x21a5c2c8ee9feb5c4b9a753b7b476b7fd6422ef1f3dd674761fa99d6ac27c8a9a197b272822f6cd57a55aa4f50ae317b13545f,
     0x15d4860d088ddb3496b0c6064756260441cde4af1771d4db01ffe5b34e59703dc255a868a1180515603aeab60794e54bb7996a7,
     0x61b1cfab6be5f32bbfa78324ed106a7636b9c5a7bd198d0158aa4f5488d08f38514f1fdf4b4f40d2181b3681c364ba0273c706,
     0x10000000000000000000000000000000000000000000000000001e2aad6a612f33307be5fa47c3c9e052f838164cd37d9a21173),
    ('sect571k1 (K-571)', '2', 4, 0x80000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000425, 0x0, 0x1,
     0x26eb7a859923fbc82189631f8103fe4ac9ca2970012d5d46024804801841ca44370958493b205e647da304db4ceb08cbbd1ba39494776fb988b47174dca88c7e2945283a01c8972,
     0x349dc807f4fbf374f4aeade3bca95314dd58cec9f307a54ffc61efc006d8a2c9d4979c0ac44aea74fbebbb9f772aedcb620b01a7ba7af1b320430c8591984f601cd4c143ef1c7a3,
     0x20000000000000000000000000000000000000000000000000000000000000000000000131850e1f19a63e4b391a8db917f4138b630d84be5d639381e91deb45cfe778f637c1001),
    ('sect571r1 (B-571)', '2', 2, 0x80000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000425, 0x1, 0x2f40e7e2221f295de297117b7f3d62f5c6a97ffcb8ceff1cd6ba8ce4a9a18ad84ffabbd8efa59332be7ad6756a66e294afd185a78ff12aa520e4de739baca0c7ffeff7f2955727a,
     0x303001d34b856296c16c0d40d3cd7750a93d1d2955fa80aa5f40fc8db7b2abdbde53950f4c0d293cdd711a35b67fb1499ae60038614f1394abfa3b4c850d927e1e7769c8eec2d19,
     0x37bf27342da639b6dccfffeb73d69d78c6c27a6009cbbca1980f8533921e8a684423e43bab08a576291af8f461bb2a8b3531d2f0485c19b16e2f1516e23dd3c1a4827af1b8ac15b,
     0x3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe661ce18ff55987308059b186823851ec7dd9ca1161de93d5174d66e8382e9bb2fe84e47),
    ('c2pnb163v1', '2', 2, 0x80000000000000000000000000000000000000107, 0x72546b5435234a422e0789675f432c89435de5242, 0xc9517d06d5240d3cff38c74b20b6cd4d6f9dd4d9,
     0x7af69989546103d79329fcc3d74880f33bbe803cb,
     0x1ec23211b5966adea1d3f87f7ea5848aef0b7ca9f,
     0x400000000000000000001e60fc8821cc74daeafc1),
    ('c2pnb163v2', '2', 2, 0x80000000000000000000000000000000000000107, 0x108b39e77c4b108bed981ed0e890e117c511cf072, 0x667aceb38af4e488c407433ffae4f1c811638df20,
     0x24266e4eb5106d0a964d92c4860e2671db9b6cc5,
     0x79f684ddf6684c5cd258b3890021b2386dfd19fc5,
     0x3fffffffffffffffffffdf64de1151adbb78f10a7),
    ('c2pnb163v3', '2', 2, 0x80000000000000000000000000000000000000107, 0x7a526c63d3e25a256a007699f5447e32ae456b50e, 0x3f7061798eb99e238fd6f1bf95b48feeb4854252b,
     0x2f9f87b7c574d0bdecf8a22e6524775f98cdebdcb,
     0x5b935590c155e17ea48eb3ff3718b893df59a05d0,
     0x3fffffffffffffffffffe1aee140f110aff961309),
    ('c2tnb191v1', '2', 2, 0x800000000000000000000000000000000000000000000201, 0x2866537b676752636a68f56554e12640276b649ef7526267, 0x2e45ef571f00786f67b0081b9495a3d95462f5de0aa185ec,
     0x36b3daf8a23206f9c4f299d7b21a9c369137f2c84ae1aa0d,
     0x765be73433b3f95e332932e70ea245ca2418ea0ef98018fb,
     0x40000000000000000000000004a20e90c39067c893bbb9a5),
    ('c2tnb191v2', '2', 4, 0x800000000000000000000000000000000000000000000201, 0x401028774d7777c7b7666d1366ea432071274f89ff01e718, 0x620048d28bcbd03b6249c99182b7c8cd19700c362c46a01,
     0x3809b2b7cc1b28cc5a87926aad83fd28789e81e2c9e3bf10,
     0x17434386626d14f3dbf01760d9213a3e1cf37aec437d668a,
     0x20000000000000000000000050508cb89f652824e06b8173),
    ('c2tnb191v3', '2', 6, 0x800000000000000000000000000000000000000000000201, 0x6c01074756099122221056911c77d77e77a777e7e7e77fcb, 0x71fe1af926cf847989efef8db459f66394d90f32ad3f15e8,
     0x375d4ce24fde434489de8746e71786015009e66e38a926dd,
     0x545a39176196575d985999366e6ad34ce0a77cd7127b06be,
     0x155555555555555555555555610c0b196812bfb6288a3ea3),
    ('c2tnb239v1', '2', 4, 0x800000000000000000000000000000000000000000000000001000000001, 0x32010857077c5431123a46b808906756f543423e8d27877578125778ac76, 0x790408f2eedaf392b012edefb3392f30f4327c0ca3f31fc383c422aa8c16,
     0x57927098fa932e7c0a96d3fd5b706ef7e5f5c156e16b7e7c86038552e91d,
     0x61d8ee5077c33fecf6f1a16b268de469c3c7744ea9a971649fc7a9616305,
     0x2000000000000000000000000000000f4d42ffe1492a4993f1cad666e447),
    ('c2tnb239v2', '2', 6, 0x800000000000000000000000000000000000000000000000001000000001, 0x4230017757a767fae42398569b746325d45313af0766266479b75654e65f, 0x5037ea654196cff0cd82b2c14a2fcf2e3ff8775285b545722f03eacdb74b,
     0x28f9d04e900069c8dc47a08534fe76d2b900b7d7ef31f5709f200c4ca205,
     0x5667334c45aff3b5a03bad9dd75e2c71a99362567d5453f7fa6e227ec833,
     0x1555555555555555555555555555553c6f2885259c31e3fcdf154624522d),
    ('c2tnb239v3', '2', 10, 0x800000000000000000000000000000000000000000000000001000000001, 0x1238774666a67766d6676f778e676b66999176666e687666d8766c66a9f, 0x6a941977ba9f6a435199acfc51067ed587f519c5ecb541b8e44111de1d40,
     0x70f6e9d04d289c4e89913ce3530bfde903977d42b146d539bf1bde4e9c92,
     0x2e5a0eaf6e5e1305b9004dce5c0ed7fe59a35608f33837c816d80b79f461,
     0xcccccccccccccccccccccccccccccac4912d2d9df903ef9888b8a0e4cff),
    ('brainpoolp160r1', 'p', 1, 0xe95e4a5f737059dc60dfc7ad95b3d8139515620f, 0x340e7be2a280eb74e2be61bada745d97e8f7c300, 0x1e589a8595423412134faa2dbdec95c8d8675e58,
     0xbed5af16ea3f6a4f62938c4631eb5af7bdbcdbc3,
     0x1667cb477a1a8ec338f94741669c976316da6321,
     0xe95e4a5f737059dc60df5991d45029409e60fc09),
    ('brainpoolp160t1', 'p', 1, 0xe95e4a5f737059dc60dfc7ad95b3d8139515620f, 0xe95e4a5f737059dc60dfc7ad95b3d8139515620c, 0x7a556b6dae535b7b51ed2c4d7daa7a0b5c55f380,
     0xb199b13b9b34efc1397e64baeb05acc265ff2378,
     0xadd6718b7c7c1961f0991b842443772152c9e0ad,
     0xe95e4a5f737059dc60df5991d45029409e60fc09),
    ('brainpoolp192r1', 'p', 1, 0xc302f41d932a36cda7a3463093d18db78fce476de1a86297, 0x6a91174076b1e0e19c39c031fe8685c1cae040e5c69a28ef, 0x469a28ef7c28cca3dc721d044f4496bcca7ef4146fbf25c9,
     0xc0a0647eaab6a48753b033c56cb0f0900a2f5c4853375fd6,
     0x14b690866abd5bb88b5f4828c1490002e6773fa2fa299b8f,
     0xc302f41d932a36cda7a3462f9e9e916b5be8f1029ac4acc1),
    ('brainpoolp192t1', 'p', 1, 0xc302f41d932a36cda7a3463093d18db78fce476de1a86297, 0xc302f41d932a36cda7a3463093d18db78fce476de1a86294, 0x13d56ffaec78681e68f9deb43b35bec2fb68542e27897b79,
     0x3ae9e58c82f63c30282e1fe7bbf43fa72c446af6f4618129,
     0x97e2c5667c2223a902ab5ca449d0084b7e5b3de7ccc01c9,
     0xc302f41d932a36cda7a3462f9e9e916b5be8f1029ac4acc1),
    ('brainpoolp224r1', 'p', 1, 0xd7c134aa264366862a18302575d1d787b09f075797da89f57ec8c0ff, 0x68a5e62ca9ce6c1c299803a6c1530b514e182ad8b0042a59cad29f43, 0x2580f63ccfe44138870713b1a92369e33e2135d266dbb372386c400b,
     0xd9029ad2c7e5cf4340823b2a87dc68c9e4ce3174c1e6efdee12c07d,
     0x58aa56f772c0726f24c6b89e4ecdac24354b9e99caa3f6d3761402cd,
     0xd7c134aa264366862a18302575d0fb98d116bc4b6ddebca3a5a7939f),
    ('brainpoolp224t1', 'p', 1, 0xd7c134aa264366862a18302575d1d787b09f075797da89f57ec8c0ff, 0xd7c134aa264366862a18302575d1d787b09f075797da89f57ec8c0fc, 0x4b337d934104cd7bef271bf60ced1ed20da14c08b3bb64f18a60888d,
     0x6ab1e344ce25ff3896424e7ffe14762ecb49f8928ac0c76029b4d580,
     0x374e9f5143e568cd23f3f4d7c0d4b1e41c8cc0d1c6abd5f1a46db4c,
     0xd7c134aa264366862a18302575d0fb98d116bc4b6ddebca3a5a7939f),
    ('brainpoolp256r1', 'p', 1, 0xa9fb57dba1eea9bc3e660a909d838d726e3bf623d52620282013481d1f6e5377, 0x7d5a0975fc2c3057eef67530417affe7fb8055c126dc5c6ce94a4b44f330b5d9, 0x26dc5c6ce94a4b44f330b5d9bbd77cbf958416295cf7e1ce6bccdc18ff8c07b6,
     0x8bd2aeb9cb7e57cb2c4b482ffc81b7afb9de27e1e3bd23c23a4453bd9ace3262,
     0x547ef835c3dac4fd97f8461a14611dc9c27745132ded8e545c1d54c72f046997,
     0xa9fb57dba1eea9bc3e660a909d838d718c397aa3b561a6f7901e0e82974856a7),
    ('brainpoolp256t1', 'p', 1, 0xa9fb57dba1eea9bc3e660a909d838d726e3bf623d52620282013481d1f6e5377, 0xa9fb57dba1eea9bc3e660a909d838d726e3bf623d52620282013481d1f6e5374, 0x662c61c430d84ea4fe66a7733d0b76b7bf93ebc4af2f49256ae58101fee92b04,
     0xa3e8eb3cc1cfe7b7732213b23a656149afa142c47aafbc2b79a191562e1305f4,
     0x2d996c823439c56d7f7b22e14644417e69bcb6de39d027001dabe8f35b25c9be,
     0xa9fb57dba1eea9bc3e660a909d838d718c397aa3b561a6f7901e0e82974856a7),
    ('brainpoolp320r1', 'p', 1, 0xd35e472036bc4fb7e13c785ed201e065f98fcfa6f6f40def4f92b9ec7893ec28fcd412b1f1b32e27, 0x3ee30b568fbab0f883ccebd46d3f3bb8a2a73513f5eb79da66190eb085ffa9f492f375a97d860eb4, 0x520883949dfdbc42d3ad198640688a6fe13f41349554b49acc31dccd884539816f5eb4ac8fb1f1a6,
     0x43bd7e9afb53d8b85289bcc48ee5bfe6f20137d10a087eb6e7871e2a10a599c710af8d0d39e20611,
     0x14fdd05545ec1cc8ab4093247f77275e0743ffed117182eaa9c77877aaac6ac7d35245d1692e8ee1,
     0xd35e472036bc4fb7e13c785ed201e065f98fcfa5b68f12a32d482ec7ee8658e98691555b44c59311),
    ('brainpoolp320t1', 'p', 1, 0xd35e472036bc4fb7e13c785ed201e065f98fcfa6f6f40def4f92b9ec7893ec28fcd412b1f1b32e27, 0xd35e472036bc4fb7e13c785ed201e065f98fcfa6f6f40def4f92b9ec7893ec28fcd412b1f1b32e24, 0xa7f561e038eb1ed560b3d147db782013064c19f27ed27c6780aaf77fb8a547ceb5b4fef422340353,
     0x925be9fb01afc6fb4d3e7d4990010f813408ab106c4f09cb7ee07868cc136fff3357f624a21bed52,
     0x63ba3a7a27483ebf6671dbef7abb30ebee084e58a0b077ad42a5a0989d1ee71b1b9bc0455fb0d2c3,
     0xd35e472036bc4fb7e13c785ed201e065f98fcfa5b68f12a32d482ec7ee8658e98691555b44c59311),
    ('brainpoolp384r1', 'p', 1, 0x8cb91e82a3386d280f5d6f7e50e641df152f7109ed5456b412b1da197fb71123acd3a729901d1a71874700133107ec53, 0x7bc382c63d8c150c3c72080ace05afa0c2bea28e4fb22787139165efba91f90f8aa5814a503ad4eb04a8c7dd22ce2826, 0x4a8c7dd22ce28268b39b55416f0447c2fb77de107dcd2a62e880ea53eeb62d57cb4390295dbc9943ab78696fa504c11,
     0x1d1c64f068cf45ffa2a63a81b7c13f6b8847a3e77ef14fe3db7fcafe0cbd10e8e826e03436d646aaef87b2e247d4af1e,
     0x8abe1d7520f9c2a45cb1eb8e95cfd55262b70b29feec5864e19c054ff99129280e4646217791811142820341263c5315,
     0x8cb91e82a3386d280f5d6f7e50e641df152f7109ed5456b31f166e6cac0425a7cf3ab6af6b7fc3103b883202e9046565),
    ('brainpoolp384t1', 'p', 1, 0x8cb91e82a3386d280f5d6f7e50e641df152f7109ed5456b412b1da197fb71123acd3a729901d1a71874700133107ec53, 0x8cb91e82a3386d280f5d6f7e50e641df152f7109ed5456b412b1da197fb71123acd3a729901d1a71874700133107ec50, 0x7f519eada7bda81bd826dba647910f8c4b9346ed8ccdc64e4b1abd11756dce1d2074aa263b88805ced70355a33b471ee,
     0x18de98b02db9a306f2afcd7235f72a819b80ab12ebd653172476fecd462aabffc4ff191b946a5f54d8d0aa2f418808cc,
     0x25ab056962d30651a114afd2755ad336747f93475b7a1fca3b88f2b6a208ccfe469408584dc2b2912675bf5b9e582928,
     0x8cb91e82a3386d280f5d6f7e50e641df152f7109ed5456b31f166e6cac0425a7cf3ab6af6b7fc3103b883202e9046565),
    ('brainpoolp512r1', 'p', 1, 0xaadd9db8dbe9c48b3fd4e6ae33c9fc07cb308db3b3c9d20ed6639cca703308717d4d9b009bc66842aecda12ae6a380e62881ff2f2d82c68528aa6056583a48f3, 0x7830a3318b603b89e2327145ac234cc594cbdd8d3df91610a83441caea9863bc2ded5d5aa8253aa10a2ef1c98b9ac8b57f1117a72bf2c7b9e7c1ac4d77fc94ca, 0x3df91610a83441caea9863bc2ded5d5aa8253aa10a2ef1c98b9ac8b57f1117a72bf2c7b9e7c1ac4d77fc94cadc083e67984050b75ebae5dd2809bd638016f723,
     0x81aee4bdd82ed9645a21322e9c4c6a9385ed9f70b5d916c1b43b62eef4d0098eff3b1f78e2d0d48d50d1687b93b97d5f7c6d5047406a5e688b352209bcb9f822,
     0x7dde385d566332ecc0eabfa9cf7822fdf209f70024a57b1aa000c55b881f8111b2dcde494a5f485e5bca4bd88a2763aed1ca2b2fa8f0540678cd1e0f3ad80892,
     0xaadd9db8dbe9c48b3fd4e6ae33c9fc07cb308db3b3c9d20ed6639cca70330870553e5c414ca92619418661197fac10471db1d381085ddaddb58796829ca90069),
    ('brainpoolp512t1', 'p', 1, 0xaadd9db8dbe9c48b3fd4e6ae33c9fc07cb308db3b3c9d20ed6639cca703308717d4d9b009bc66842aecda12ae6a380e62881ff2f2d82c68528aa6056583a48f3, 0xaadd9db8dbe9c48b3fd4e6ae33c9fc07cb308db3b3c9d20ed6639cca703308717d4d9b009bc66842aecda12ae6a380e62881ff2f2d82c68528aa6056583a48f0, 0x7cbbbcf9441cfab76e1890e46884eae321f70c0bcb4981527897504bec3e36a62bcdfa2304976540f6450085f2dae145c22553b465763689180ea2571867423e,
     0x640ece5c12788717b9c1ba06cbc2a6feba85842458c56dde9db1758d39c0313d82ba51735cdb3ea499aa77a7d6943a64f7a3f25fe26f06b51baa2696fa9035da,
     0x5b534bd595f5af0fa2c892376c84ace1bb4e3019b71634c01131159cae03cee9d9932184beef216bd71df2dadf86a627306ecff96dbb8bace198b61e00f8b332,
     0xaadd9db8dbe9c48b3fd4e6ae33c9fc07cb308db3b3c9d20ed6639cca70330870553e5c414ca92619418661197fac10471db1d381085ddaddb58796829ca90069),
    ('sm2p256v1', 'p', 1, 0xfffffffeffffffffffffffffffffffffffffffff00000000ffffffffffffffff, 0xfffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffc, 0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93,
     0x32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7,
     0xbc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0,
     0xfffffffeffffffffffffffffffffffff7203df6b21c6052b53bbf40939d54123),
    ('c2tnb359v1', '2', 76, 0x800000000000000000000000000000000000000000000000000000000000000000000000100000000000000001, 0x5667676a654b20754f356ea92017d946567c46675556f19556a04616b567d223a5e05656fb549016a96656a557, 0x2472e2d0197c49363f1fe7f5b6db075d52b6947d135d8ca445805d39bc345626089687742b6329e70680231988,
     0x3c258ef3047767e7ede0f1fdaa79daee3841366a132e163aced4ed2401df9c6bdcde98e8e707c07a2239b1b097,
     0x53d7e08529547048121e9c95f3791dd804963948f34fae7bf44ea82365dc7868fe57e4ae2de211305a407104bd,
     0x1af286bca1af286bca1af286bca1af286bca1af286bc9fb8f6b85c556892c20a7eb964fe7719e74f490758d3b),
]

# checked only directly (cofactor 10080 is too large for the blind tests here)
LARGE = [
    ('c2tnb431r1', '2', 10080, 0x800000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000001, 0x1a827ef00dd6fc0e234caf046c6a5d8a85395b236cc4ad2cf32a0cadbdc9ddf620b0eb9906d0957f6c6feacd615468df104de296cd8f, 0x10d9b4a3d9047d8b154359abfb1b7f5485b04ceb868237ddc9deda982a679a5a919b626d4e50a8dd731b107a9962381fb5d807bf2618,
     0x120fc05d3c67a99de161d2f4092622feca701be4f50f4758714e8a87bbf2a658ef8c21e7c5efe965361f6c2999c0c247b0dbd70ce6b7,
     0x20d0af8903a96f8d5fa2c255745d3c451b302c9346d9b7e485e7bce41f6b591f3e8f6addcbb0bc4c2f947a7de1a89b625d6a598b3760,
     0x340340340340340340340340340340340340340340340340340340323c313fab50589703b5ec68d3587fec60d161cc149c1ad4a91),
]

# ---------------------------------------------------------------- arithmetic
def gf2_mul(a, b, f, m):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if (a >> m) & 1: a ^= f
    return r

def gf2_inv(a, f):
    u, v, g1, g2 = a, f, 1, 0
    while u != 1:
        j = u.bit_length() - v.bit_length()
        if j < 0: u, v, g1, g2, j = v, u, g2, g1, -j
        u ^= v << j; g1 ^= g2 << j
    return g1

class Binary:
    """y^2 + xy = x^3 + a x^2 + b over GF(2^m) = GF(2)[z]/f(z). All curves used here have odd m."""
    def __init__(s, f, a, b): s.f, s.m, s.a, s.b = f, f.bit_length() - 1, a, b
    def M(s, x, y): return gf2_mul(x, y, s.f, s.m)
    def on(s, P): x, y = P; return s.M(y, y) ^ s.M(x, y) == s.M(s.M(x, x), x) ^ s.M(s.a, s.M(x, x)) ^ s.b
    def add(s, P, Q):
        if P is None: return Q
        if Q is None: return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if y1 ^ y2 == x1 or x1 == 0: return None
            l = x1 ^ s.M(y1, gf2_inv(x1, s.f)); x3 = s.M(l, l) ^ l ^ s.a
            return (x3, s.M(x1, x1) ^ s.M(l ^ 1, x3))
        l = s.M(y1 ^ y2, gf2_inv(x1 ^ x2, s.f)); x3 = s.M(l, l) ^ l ^ x1 ^ x2 ^ s.a
        return (x3, s.M(l, x1 ^ x3) ^ x3 ^ y1)
    def lift(s, x):
        """A point with this x, or None. Solves z^2 + z = (x^3 + a x^2 + b)/x^2 with the half-trace (odd m)."""
        if x == 0: return None
        v = s.M(s.M(s.M(x, x), x) ^ s.M(s.a, s.M(x, x)) ^ s.b, gf2_inv(s.M(x, x), s.f))
        t, w = v, v
        for _ in range(s.m - 1): w = s.M(w, w); t ^= w
        if t: return None                      # trace 1: no point with this x
        z, w = v, v
        for _ in range((s.m - 1) // 2): w = s.M(s.M(w, w), s.M(w, w)); z ^= w
        return (x, s.M(x, z))

def sqrt_mod(v, p):
    v %= p
    if v == 0: return 0
    if pow(v, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(v, (p + 1) // 4, p)
    q, e = p - 1, 0
    while q % 2 == 0: q //= 2; e += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = e, pow(z, q, p), pow(v, q, p), pow(v, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

class Prime:
    """y^2 = x^3 + a x + b over GF(p)."""
    def __init__(s, p, a, b): s.p, s.a, s.b, s.m = p, a, b, p.bit_length()
    def on(s, P): x, y = P; return (y * y - x ** 3 - s.a * x - s.b) % s.p == 0
    def add(s, P, Q):
        if P is None: return Q
        if Q is None: return P
        p = s.p; (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0: return None
            l = (3 * x1 * x1 + s.a) * pow(2 * y1, -1, p) % p
        else:
            l = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (l * l - x1 - x2) % p
        return (x3, (l * (x1 - x3) - y1) % p)
    def lift(s, x):
        y = sqrt_mod(x ** 3 + s.a * x + s.b, s.p)
        return None if y is None else (x, y)

def mul(E, k, P):
    R = None
    while k:
        if k & 1: R = E.add(R, P)
        P = E.add(P, P); k >>= 1
    return R

def primes_of(h):
    out, d = [], 2
    while d * d <= h:
        if h % d == 0: out.append(d)
        while h % d == 0: h //= d
        d += 1
    return out + ([h] if h > 1 else [])

def torsion_generator(E, n, h):
    """A point of order exactly h (deterministic: lift x = 1, 2, 3, ... and multiply by n). All cofactor groups here are cyclic."""
    x = 0
    while True:
        x += 1
        P = E.lift(x)
        if P is None: continue
        T = mul(E, n, P)
        assert mul(E, h, T) is None, 'cofactor mismatch'
        if all(mul(E, h // q, T) is not None for q in primes_of(h)): return T

def hidden_points(E, G, n, h):
    """All (k, P) with k*P = G, for every divisor k > 1 of h (k = 2 when h = 1)."""
    if h == 1: return [(2, mul(E, pow(2, -1, n), G))]
    T = torsion_generator(E, n, h); out = []
    for k in [d for d in range(2, h + 1) if h % d == 0]:
        P, step = mul(E, pow(k, -1, n), G), mul(E, h // k, T)
        for _ in range(k):                      # P + j*(h/k)*T runs over every solution of k*P = G
            if P is not None:
                assert mul(E, k, P) == G
                out.append((k, P))
            P = E.add(P, step)
    return out

def longest_run(a, b):
    """Length of the longest common substring of two strings."""
    best, prev = 0, [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best: best = cur[j]
        prev = cur
    return best

def chi2(s):
    """Chi-square of the 16 hex-digit counts against uniform (15 degrees of freedom)."""
    e = len(s) / 16
    return sum((s.count(c) - e) ** 2 / e for c in '0123456789abcdef')

# ---------------------------------------------------------------- run
def main():
    H, PUB = {}, {}
    for name, t, h, p, a, b, gx, gy, n in CURVES:
        E = Prime(p, a, b) if t == 'p' else Binary(p, a, b)
        G = (gx, gy)
        assert E.on(G) and mul(E, n, G) is None, name + ': generator check failed'
        H[name] = (E.m, [(k, P[0]) for k, P in hidden_points(E, G, n, h)])
        PUB[name] = gx
    total = sum(len(v[1]) for v in H.values())
    print(f'{len(CURVES)} curves, every generator verified on-curve with order n; {total} hidden points computed.\n')

    print('TEST 1  Hidden x-coordinates much shorter than chance (shortfall >= 20 bits below field size - log2(count))')
    for name, (m, xs) in H.items():
        x = min(v for _, v in xs); short = m - math.log2(len(xs)) - x.bit_length()
        if short >= 20:
            k = next(k for k, v in xs if v == x)
            print(f'   {name:18} G = {k}P   x(P) = 0x{x:x}   ({x.bit_length()} bits on a {m}-bit field)')
    print('   (all except sect113r2 were reported before; see the note)\n')

    print('TEST 2  Blind, every pair of curves: longest run of hex digits shared by their hidden x-coordinates')
    hist, hits = {}, []
    for (a, (_, xa)), (b, (_, xb)) in combinations(H.items(), 2):
        L, ka, sa, kb, sb = max((longest_run(format(u, 'x'), format(v, 'x')), k1, u, k2, v) for k1, u in xa for k2, v in xb)
        hist[L] = hist.get(L, 0) + 1
        if L >= 10: hits.append((a, ka, sa, b, kb, sb, L))
    for a, ka, sa, b, kb, sb, L in hits:
        print(f'   {L:2} digits  {a} (G = {ka}P)  and  {b} (G = {kb}P)\n              0x{sa:x}\n              0x{sb:x}')
    print('   pairs by longest shared run:', ', '.join(f'{k}: {v}' for k, v in sorted(hist.items())), '\n')

    print("TEST 3  Blind: hidden x-coordinates that share a long run with another curve's PUBLISHED generator x")
    hist3 = {}
    for name, (_, xs) in H.items():
        for k, v in xs:
            for other, gx in PUB.items():
                if other == name: continue
                L = longest_run(format(v, 'x'), format(gx, 'x')); hist3[L] = hist3.get(L, 0) + 1
                if L >= 10:
                    print(f'   {L:2} digits  hidden x of {name} (G = {k}P)  vs  published Gx of {other}\n              0x{v:x}\n              0x{gx:x}')
    print('   comparisons by longest shared run:', ', '.join(f'{k}: {v}' for k, v in sorted(hist3.items())), '\n')

    print('TEST 4  Blind: hex-digit unevenness (chi-square, 15 d.f.); random strings of these lengths essentially never exceed 60')
    import random
    rng = random.Random(1)
    rand_max = max(chi2(''.join(rng.choice('0123456789abcdef') for _ in range(48))) for _ in range(20000))
    flagged = []
    for name, (_, xs) in H.items():
        for k, v in xs:
            if chi2(format(v, 'x')) > 60: flagged.append((chi2(format(v, 'x')), f'hidden x of {name} (G = {k}P)', v))
        if chi2(format(PUB[name], 'x')) > 50: flagged.append((chi2(format(PUB[name], 'x')), f'published Gx of {name}', PUB[name]))
    for c2, what, v in sorted(flagged, reverse=True):
        print(f'   chi2 {c2:5.1f}  {what}\n              0x{v:x}')
    print(f'   (largest chi-square among 20,000 random 48-digit strings: {rand_max:.1f})\n')

    # ------------------------------------------------------------ the note's specific claims
    def direct_check(row, k, x):
        name, t, h, p, a, b, gx, gy, n = row
        E = Binary(p, a, b); G = (gx, gy); P = E.lift(x)
        return P is not None and mul(E, k, P) in (G, (gx, gx ^ gy))     # -P maps to -G = (x, x + y)
    X431 = 0x38f789767b657576e655447567670102792657bec8975851349579234562399235492978123847a777c67666f6f6766a666666b66610
    def hx(name, k): return {format(v, 'x') for kk, v in H[name][1] if kk == k}
    B = '388477117474774774375771981418781cabe047772771'
    claims = [
        ('sect113r1: G = 2P with x(P) = 3',                     '3' in hx('sect113r1', 2)),
        ('sect113r2: G = 2P with x(P) = 2',                     '2' in hx('sect113r2', 2)),
        ('sect163k1 (NIST K-163): 32-hex-digit x(P), G = 2P',         '23e21d6019e1211f6bd47ec180256e97' in hx('sect163k1 (K-163)', 2)),
        ('c2pnb163v1: 32-hex-digit x(P), G = 2P',                     '13152a6b92ceb34c0969c282324645cb' in hx('c2pnb163v1', 2)),
        ('sect131r1 / r2: one 31-digit body, last digit 2 / 4',  '10134db223421eb623422888a7879682' in hx('sect131r1', 2) and '10134db223421eb623422888a7879684' in hx('sect131r2', 2)),
        ('sect193r1 and sect193r2: identical x(P)',              hx('sect193r1', 2) & hx('sect193r2', 2) == {'23daf8a23206f9c4f299d7b21a9c369137f2c84ae1aa04'}),
        ('c2tnb191v1 G = 2P, x(P) = 3.B.4',                      '3' + B + '4' in hx('c2tnb191v1', 2)),
        ('c2tnb191v2 G = 4P, x(P) = 3.B.d',                      '3' + B + 'd' in hx('c2tnb191v2', 4)),
        ('c2tnb191v3 G = 6P, x(P) = 3.B.0',                      '3' + B + '0' in hx('c2tnb191v3', 6)),
        ('c2tnb239v1 G = 4P, x(P) = 5757157563.B.2',             '5757157563' + B + '2' in hx('c2tnb239v1', 4)),
        ('c2tnb239v2 G = 6P, x(P) = 5757157563.B.4',             '5757157563' + B + '4' in hx('c2tnb239v2', 6)),
        ('secp112r2 G = 4P, x(P) = published Gx of secp112r1',   format(PUB['secp112r1'], 'x') in hx('secp112r2', 4)),
        ('secp128r2 G = 4P, x(P) = published Gx of secp128r1',   format(PUB['secp128r1'], 'x') in hx('secp128r2', 4)),
        ('c2pnb163v3 G = 2P, x(P) = published Gx of c2pnb163v1', format(PUB['c2pnb163v1'], 'x') in hx('c2pnb163v3', 2)),
        ('c2pnb163v2 G = 2P, x(P) = c2pnb163v1 Gx less one digit', '7f69989546103d79329fcc3d74880f33bbe803cb' in hx('c2pnb163v2', 2)),
        ('secp160r2 G = 2P, x(P) starts with 23 digits of secp160r1 Gx', any(v.startswith(format(PUB['secp160r1'], 'x')[:23]) for v in hx('secp160r2', 2))),
        ('sect193r1/r2 x(P) contains 44 digits of c2tnb191v1 Gx', format(PUB['c2tnb191v1'], 'x')[3:-1] in next(iter(hx('sect193r1', 2) & hx('sect193r2', 2)))),
        ('c2tnb359v1 G = 76P, 129-bit x(P) sharing 29 digits with c2pnb163v1\'s hidden x',
                                                                 '152d52a6b92ceb34c0969c282324645cb' in hx('c2tnb359v1', 76)),
        ('c2tnb239v3 G = 10P, x(P) with highly uneven digits (chi-square > 90)',
                                                                 any(chi2(v) > 90 for v in hx('c2tnb239v3', 10))),
        ('c2tnb431r1 G = 10080P for the patterned x(P) 38f78976...b66610', direct_check(LARGE[0], 10080, X431)),
        ('secp256k1 (Bitcoin) G = 2P, x(P) = 0x3b78ce...9c63',   '3b78ce563f89a0ed9414f5aa28ad0d96d6795f9c63' in hx('secp256k1', 2)),
    ]
    print('CLAIMS IN THE NOTE')
    for text, ok in claims: print(f'   {"PASS" if ok else "FAIL"}  {text}')
    bad = sum(not ok for _, ok in claims)
    print(f'\n{len(claims) - bad} of {len(claims)} claims pass.')
    raise SystemExit(1 if bad else 0)

if __name__ == '__main__':
    main()
