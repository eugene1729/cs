import Data.List
import Data.Ratio
import qualified Data.Vector as V

[] .+. ys = ys
xs .+. [] = xs
(x:xs) .+. (y:ys) = (x + y) : xs .+. ys
infixl 6 .+.

[] .*. _ = []
_ .*. [] = []
(0:p1) .*. p2 = 0 : p1 .*. p2
p1 .*. (0:p2) = 0 : p1 .*. p2
(1:p1) .*. p2 = p2 .+. (0 : p1 .*. p2)
p1 .*. (1:p2) = p1 .+. (0 : p1 .*. p2)
(p:p1) .*. p2 = map (p *) p2 .+. (0 : p1 .*. p2)
infixl 7 .*.

c n k | k == 0 = 1
      | k == n = 1
      | k > n = 0
      | otherwise = foldl (\a i -> (a * (n - i + 1)) `div` i) 1 [1 .. k]

fac = product . enumFromTo 1

sortDescending :: Ord a => [a] -> [a]
sortDescending = sortBy (flip compare)

change :: [Int] -> Integer -> Integer
change coins amount | amount > 0 && amount `mod` g' == 0 = sum (map f [0, l .. t])
                    | otherwise = 0 where
    g' = fromIntegral g
    g = foldr1 gcd coins
    cs = map (`div` g) (sortDescending coins)
    a = amount `div` g'
    l = foldr1 lcm cs
    l' = fromIntegral l
    n = genericLength cs
    p = V.fromList $ foldr1 (.*.) $ map (\c -> take (l - c + 1) (cycle (1 : replicate (c - 1) 0))) cs
    t = V.length p
    f d | fromIntegral d <= a && i < t = c (k + n - 1) (n - 1) * p V.! i
        | otherwise = 0 where 
        k = (a - fromIntegral d) `div` l'
        i = fromIntegral (a - k * l')

poly :: [Int] -> Int -> [Rational]
poly coins r = zipWith (%) (foldr1 (.+.) $ map f [0, l .. t]) [ l' ^ i * fac (fromIntegral n - 1) | i <- [0 .. n - 1] ] where
    cs = sortDescending coins
    a = r `mod` l
    l = foldr1 lcm cs
    l' = fromIntegral l
    n = length cs
    p = V.fromList $ foldr1 (.*.) $ map (\c -> take (l - c + 1) (cycle (1 : replicate (c - 1) 0))) cs
    t = V.length p
    f d | i < t = map (fromIntegral (p V.! i) *) ps
        | otherwise = [] where 
        i = a - ((a - d) `div` l) * l
        ps | n >= 2 = foldr1 (.*.) (map (\j -> [ fromIntegral (n - j - d `div` l), 1]) [1 .. n - 1])
           | otherwise = [1]

usCoins = [1, 5, 10, 25, 50, 100] :: [Int] -- in cents

usNotes = [1, 2, 5, 10, 20, 50, 100] :: [Int] -- in dollars

usCoinsAndNotes = usCoins ++ (map (100 *) usNotes)

googol = 100 * 10 ^ 100 -- in cents
