function tests = tests()
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.TestData = dome([],[],[],'localfunctions');
end

function test_rotmat(testCase)
    rotmat = testCase.TestData.rotmat;
    verifyEqual(testCase, rotmat(90)*[1;0], [0;1])
    verifyEqual(testCase, rotmat(-90)*[1;1], [1;-1])
end

function test_reflection_matrix(testCase)
    reflection_matrix = testCase.TestData.reflection_matrix;
    verifyEqual(testCase, reflection_matrix(0, 1)*[1;0], [-1;0])
    verifyEqual(testCase, reflection_matrix(1, 1)*[1;0], [0;1])
end

function test_reflect(testCase)
    reflect = testCase.TestData.reflect;
    verifyEqual(testCase, reflect([1,0], [0,0], [0,1]), [-1, 0])
    verifyEqual(testCase, reflect([1,0; 0,2], [1,1], [2,2]), [0,1; 2,0])
end

function test_polygon_arc(testCase)
    polygon_arc = testCase.TestData.polygon_arc;
    verifyEqual(testCase, polygon_arc([-45, 45], [-90, 0, 90]), [0, sqrt(2), 2*sqrt(2)])
    verifyEqual(testCase, polygon_arc([-45, 45], [-90, 0, 90]'), [0, sqrt(2), 2*sqrt(2)]')
end

function test_wedge_faces(testCase)
    wedge_faces = testCase.TestData.wedge_faces;
    verifyEqual(testCase, wedge_faces(1, 6, 4), [6 1 3 NaN;3 4 7 6;4 5 8 7;5 2 8 NaN])
    verifyEqual(testCase, wedge_faces(4, 6, 4), [3 1 12 NaN;12 13 4 3;13 14 5 4;14 2 5 NaN])
end

function test_flap(testCase)
    flap = testCase.TestData.flap;
    [p, q, r, s] = deal([0,0], [1,0], [1,1], [0,1]);

    c = [-1,0];
    verifyEqual(testCase, flap(p, q, r, s, c, 0), [0,-1; 1,-1]);

    c = [-2,1];
    verifyEqual(testCase, flap(p, q, r, s, c, 0), [0,-1; 0.25,-1], 'abstol', 2*eps);
    
    c = [3,1];
    verifyEqual(testCase, flap(p, q, r, s, c, 0), [0.75,-1; 1,-1], 'abstol', 2*eps);
end

function test_intersect(testCase)
    intersect = testCase.TestData.intersect;
    verifyEqual(testCase, intersect([0,0], [1,1], [1,0], [0,1]), [0.5,0.5]);
end