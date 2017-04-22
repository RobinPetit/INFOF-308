qq[n, k] :=
    block(
        [res, res1],
        if k < n-1 or k > n*(n-1)/2 then
            0
        else
            if k = n-1 then
                n^(n-2)
            else
            (
                res: binomial(n*(n-1)/2, k),
                for m: 0 step 1 thru n-2 do
                (
                    res1: 0,
                    for p: max(0, k - m*(m+1)/2) step 1 thru k-m do
                        res1: res1 + binomial((n-1-m)*(n-2-m)/2, p) * qq[m+1, k-p],
                    res: res - binomial(n-1, m) * res1
                ),
                res
            )
    );

HP[p] :=
    block(
        [res: 0],
        if p < 2 then
            1
        else
        (
            for k: 1 step 1 thru p-1 do
                res: res + binomial(p-2, k-1)*(2^k-1)*HP[k]*HP[p-k],
            res
        )
    );

for n: 2 step 1 thru 50 do
(
    res: 0,
    for k: n-1 step 1 thru n*(n-1)/2 do
        res: res + qq[n, k],
    if HP[n] # res then
      print(n)
    else
      print("ok for n ==", n)
);
