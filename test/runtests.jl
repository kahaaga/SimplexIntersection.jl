using SimplexIntersection
using Base.Test

approxtol = 1/10^8

X1 =  [[0.;0.;0.] [1.;0.;0.] [0.;1.;0.] [0.;0.;1.]]
Y1 = [0.864401 0.128064 -0.251176 -0.902914;
        2.21188 1.85278 0.369714 0.864401;
        0.532813 -0.827763 0.0721164 2.21188]


#Test 2
X2 = [0.934702 0.297735 0.670881 0.85225;
        0.670881 0.626323 0.521743 0.942819;
        0.521743 0.385638 0.202466 0.613331];
Y2 = [0.670881 0.85225 0.297735 0.202466;
        0.521743 0.942819 0.626323 0.898955;
        0.202466 0.613331 0.385638 0.85225];

# Test 4
X4 =  [[0.;0.;0.] [1.;0.;0.] [0.;1.;0.] [0.;0.;1.]];
Y4 = [[0.3;0.3;0.1] [0.4;-0.1;-0.5] [0.5;0.3;1.5] [-1;0.1;0.3]];



# Test 5
X5 =  [[0.;0.;0.] [1.;0.;0.] [0.;1.;0.] [0.;0.;1.]];
Y5 = [[1;0.5;0] [-0.2;0.5;0] [0.5;-0.2;0] [1;0.5;1.5]];


# Test 6
X6 = [0.864401 0.128064 -0.251176 -0.902914;
        2.21188 1.85278 0.369714 0.864401;
        0.532813 -0.827763 0.0721164 2.21188];
Y6 = [-0.0192918 0.128064 -0.902914 0.864401;
        0.128064 1.85278 0.864401 2.21188;
        1.85278 -0.827763 2.21188 0.53];

# Test 7
X7 = [-0.271735 0.128064 0.369714 2.21188;
        0.502334 1.85278 0.0721164 0.532813;
        -0.516984 -0.827763 -1.50343 -0.271735]
Y7 = [0.128064 2.21188 0.369714 0.864401;
        1.85278 0.532813 0.0721164 2.21188;
        -0.827763 -0.271735 -1.50343 0.532813]

# Test 8
X8 = [0.227011 0.659793 0.628922 0.939909;
        0.659793 0.843993 0.877505 0.212366;
        0.843993 0.495075 0.404161 0.628922];
Y8 = [0.659793 0.939909 0.212366 0.628922;
        0.843993 0.212366 0.628922 0.877505;
        0.495075 0.628922 0.877505 0.404161];

# Test 9
X9 = [0.934702 0.297735 0.670881 0.85225;
        0.670881 0.626323 0.521743 0.942819;
        0.521743 0.385638 0.202466 0.613331];
Y9 = [0.670881 0.85225 0.297735 0.202466;
        0.521743 0.942819 0.626323 0.898955;
        0.202466 0.613331 0.385638 0.85225];




# Test 10
X10 = [0.385638 0.297735 0.353544 0.202466;
        0.353544 0.626323 0.444382 0.898955;
        0.444382 0.385638 0.235286 0.85225];
Y10 = [0.353544 0.235286 0.297735 0.202466;
        0.444382 0.822361 0.626323 0.898955;
        0.235286 0.934702 0.385638 0.85225];

X11 = [1.12545 0.0978862 0.401808 1.91521;
        0.979468 -0.171237 0.0978862 1.12545;
        0.491734 -0.452719 -0.171237 0.979468];
Y11 = [0.979468 -0.171237 0.0978862 1.12545;
        0.491734 -0.452719 -0.171237 0.979468;
        0.401808 -0.537488 -0.452719 0.491734];

X12 =  [1.0 1.4 3.0 -3.0;
        2.0 -3.0 7.0 -5.0;
        3.0 -5.0 1.4 7.0]
Y12 = [2.0 -3.0 7.0 -5.0;
        3.0 -5.0 1.4 7.0;
        7.0 7.0 -3.0 -1.0]


X19 = [0.202338 0.311773 0.18092 0.117;
        0.458579 0.446821 0.43954 0.117;
        0.144628 0.18092 0.202338 0.341191]
Y19 =   [0.338291 0.341191 0.117 0.288399;
        0.416521 0.412559 0.341191 0.464648;
        0.264743 0.274673 0.412559 0.125069]

using SimplexIntersection
@testset "simplexintersection: compare with matlab" begin
        @test SimplexIntersection.simplexintersection(X1, Y1, what = "volume") ≈ 0.00059113 atol = approxtol
        @test SimplexIntersection.simplexintersection(X4, Y4, what = "volume") ≈ 0.271468970913881 atol = approxtol
        @test SimplexIntersection.simplexintersection(X5, Y5, what = "volume") ≈ 0.371760204081633 atol = approxtol
        @test SimplexIntersection.simplexintersection(X6, Y6, what = "volume") ≈ 2.346911657841576 atol = approxtol
        @test SimplexIntersection.simplexintersection(X7, Y7, what = "volume") ≈  1.224294204983859 atol = approxtol
        @test SimplexIntersection.simplexintersection(X11, Y11, what = "volume") ≈ 0 atol = approxtol

        @test SimplexIntersection.simplexintersection(X10, Y10, what = "volume") ≈ 0.0002709969231403119 atol = approxtol
        @test SimplexIntersection.simplexintersection(X19, Y19, what = "volume") ≈ 0.000000283415912320120 atol = approxtol
        @test SimplexIntersection.simplexintersection(X2, Y2, what = "volume") ≈ 0.010081146644461 atol = approxtol
        @test SimplexIntersection.simplexintersection(X8, Y8, what = "volume") ≈ 0.029062035528103 atol = approxtol
        @test SimplexIntersection.simplexintersection(X9, Y9, what = "volume") ≈ 0.010081146644461 atol = approxtol
end;
