include("../Polyhedron.jl") 
 
 coordinates=[ [ [ 0., 1., 0. ], [ 0., 0., 0. ], [ -1., 0., 0. ], [ -1.500000000000000222, -0.86602540378443859659, 0. ], [ -2.0000000000000004441, -1.7320508075688771932, 0. ], 
 [ -1.5000000000000008882, -2.5980762113533160118, 0. ], [ -0.63397459621556206955, -2.0980762113533164559, 0. ], [ 0.23205080756887608295, -2.5980762113533173441, 0. ], 
 [ 0.73205080756887697113, -1.7320508075688791916, 0. ], [ 1.7320508075688767491, -1.7320508075688807459, 0. ], [ 1.7320508075688783034, -0.73205080756888096793, 0. ], 
 [ 0.86602540378444070601, -0.23205080756887952464, 0. ], [ 1.3660254037844419273, 0.63397459621555818376, 0. ], [ 0.86602540378444348157, 1.4999999999999975575, 0. ], [ -0.25, 0.75, 1. ], 
 [ -0.75, -0.25000000000000011102, 1. ], [ -1.3415063509461098157, -0.091506350946109538125, 1. ], [ -1.6584936490538906284, -1.6405444566227675995, 1. ], 
 [ -2.0915063509461102598, -2.0735571585149870089, 1. ], [ -0.97548094716167188523, -2.0065698604072066402, 1. ], [ -0.54246824526945258693, -2.4395825622994262716, 1. ], 
 [ 0.39054445662276709994, -1.8235571585149887852, 1. ], [ 0.98205080756887686011, -1.9820508075688794136, 1. ], [ 1.4820508075688780814, -0.98205080756888085691, 1. ], 
 [ 1.6405444566227689318, -0.39054445662277115225, 1. ], [ 1.0245190528383321116, 0.54246824526944870115, 1. ], [ 1.457531754730551965, 0.97548094716166788842, 1. ], 
 [ 0.34150635094611381248, 0.9084936490538890741, 1. ], [ 0., 1., 1. ], [ 0., 0., 1. ], [ -1., 0., 1. ], [ -1.500000000000000222, -0.86602540378443859659, 1. ], 
 [ -2.0000000000000004441, -1.7320508075688771932, 1. ], [ -1.5000000000000008882, -2.5980762113533160118, 1. ], [ -0.63397459621556206955, -2.0980762113533164559, 1. ], 
 [ 0.23205080756887608295, -2.5980762113533173441, 1. ], [ 0.73205080756887697113, -1.7320508075688791916, 1. ], [ 1.7320508075688767491, -1.7320508075688807459, 1. ], 
 [ 1.7320508075688783034, -0.73205080756888096793, 1. ], [ 0.86602540378444070601, -0.23205080756887952464, 1. ], [ 1.3660254037844419273, 0.63397459621555818376, 1. ], 
 [ 0.86602540378444348157, 1.4999999999999975575, 1. ] ], [ [ 0.73205080756887697113, -1.7320508075688791916, 0. ], [ 0.23205080756887608295, -2.5980762113533173441, 0. ], 
 [ -0.63397459621556206955, -2.0980762113533164559, 0. ], [ -1.5000000000000015543, -2.5980762113533155677, 0. ], [ -2.3660254037844410391, -3.0980762113533142355, 0. ], 
 [ -2.3660254037844419273, -4.0980762113533142355, 0. ], [ -1.3660254037844421493, -4.0980762113533160118, 0. ], [ -0.86602540378444392566, -4.9641016151377552745, 0. ], 
 [ -4.4408920985006261617e-15, -4.4641016151377579391, 0. ], [ 0.86602540378443282343, -4.9641016151377597154, 0. ], [ 1.3660254037844348218, -4.0980762113533222291, 0. ], 
 [ 0.86602540378443748637, -3.2320508075688820782, 0. ], [ 1.7320508075688771932, -2.7320508075688847427, 0. ], [ 1.7320508075688800798, -1.7320508075688847427, 0. ], 
 [ 0.39054445662276721096, -1.8235571585149885632, 1. ], [ -0.54246824526945269795, -2.4395825622994262716, 1. ], [ -0.97548094716167166318, -2.0065698604072061961, 1. ], 
 [ -2.0245190528383307793, -3.1895825622994244952, 1. ], [ -2.6160254037844410391, -3.3480762113533142355, 1. ], [ -1.6160254037844421493, -3.8480762113533151236, 1. ], 
 [ -1.4575317547305524091, -4.4395825622994253834, 1. ], [ -0.34150635094611403453, -4.3725952641916467911, 1. ], [ 0.091506350946104708655, -4.8056079660838673107, 1. ], 
 [ 1.0245190528383250061, -4.1895825622994316006, 1. ], [ 1.4575317547305453036, -3.7565698604072128575, 1. ], [ 1.3905444566227675995, -2.6405444566227744829, 1. ], 
 [ 1.9820508075688778593, -2.4820508075688851868, 1. ], [ 0.98205080756887963567, -1.9820508075688825222, 1. ], [ 0.73205080756887697113, -1.7320508075688791916, 1. ], 
 [ 0.23205080756887608295, -2.5980762113533173441, 1. ], [ -0.63397459621556206955, -2.0980762113533164559, 1. ], [ -1.5000000000000015543, -2.5980762113533155677, 1. ], 
 [ -2.3660254037844410391, -3.0980762113533142355, 1. ], [ -2.3660254037844419273, -4.0980762113533142355, 1. ], [ -1.3660254037844421493, -4.0980762113533160118, 1. ], 
 [ -0.86602540378444392566, -4.9641016151377552745, 1. ], [ -4.4408920985006261617e-15, -4.4641016151377579391, 1. ], [ 0.86602540378443282343, -4.9641016151377597154, 1. ], 
 [ 1.3660254037844348218, -4.0980762113533222291, 1. ], [ 0.86602540378443748637, -3.2320508075688820782, 1. ], [ 1.7320508075688771932, -2.7320508075688847427, 1. ], 
 [ 1.7320508075688800798, -1.7320508075688847427, 1. ] ], [ [ 3.7320508075688731964, -4.7320508075688900718, 0. ], [ 3.2320508075688749727, -3.866025403784450365, 0. ], 
 [ 4.0980762113533142355, -3.3660254037844516972, 0. ], [ 4.0980762113533160118, -2.3660254037844516972, 0. ], [ 4.0980762113533177882, -1.3660254037844512531, 0. ], 
 [ 3.23205080756888119, -0.86602540378444947677, 0. ], [ 2.7320508075688794136, -1.7320508075688871852, 0. ], [ 1.7320508075688800798, -1.7320508075688847427, 0. ], 
 [ 1.7320508075688771932, -2.7320508075688847427, 0. ], [ 0.86602540378443748637, -3.2320508075688820782, 0. ], [ 1.3660254037844343777, -4.0980762113533222291, 0. ], 
 [ 2.3660254037844339337, -4.0980762113533257818, 0. ], [ 2.3660254037844303809, -5.0980762113533248936, 0. ], [ 3.2320508075688669791, -5.5980762113533284463, 0. ], 
 [ 3.8235571585149834561, -4.390544456622779812, 1. ], [ 3.7565698604072048639, -3.2745190528383414375, 1. ], [ 4.3480762113533151236, -3.1160254037844516972, 1. ], 
 [ 3.8480762113533177882, -1.616025403784450809, 1. ], [ 4.0065698604072093048, -1.0245190528383409934, 1. ], [ 3.0735571585149896734, -1.6405444566227780356, 1. ], 
 [ 2.4820508075688803018, -1.4820508075688865191, 1. ], [ 1.9820508075688778593, -2.4820508075688851868, 1. ], [ 1.3905444566227678216, -2.6405444566227744829, 1. ], 
 [ 1.4575317547305446375, -3.7565698604072128575, 1. ], [ 1.6160254037844334896, -4.3480762113533231172, 1. ], [ 2.616025403784430825, -4.8480762113533257818, 1. ], 
 [ 2.4575317547305388644, -5.4395825622994351534, 1. ], [ 3.3905444566227593839, -4.8235571585150003315, 1. ], [ 3.7320508075688731964, -4.7320508075688900718, 1. ], 
 [ 3.2320508075688749727, -3.866025403784450365, 1. ], [ 4.0980762113533142355, -3.3660254037844516972, 1. ], [ 4.0980762113533160118, -2.3660254037844516972, 1. ], 
 [ 4.0980762113533177882, -1.3660254037844512531, 1. ], [ 3.23205080756888119, -0.86602540378444947677, 1. ], [ 2.7320508075688794136, -1.7320508075688871852, 1. ], 
 [ 1.7320508075688800798, -1.7320508075688847427, 1. ], [ 1.7320508075688771932, -2.7320508075688847427, 1. ], [ 0.86602540378443748637, -3.2320508075688820782, 1. ], 
 [ 1.3660254037844343777, -4.0980762113533222291, 1. ], [ 2.3660254037844339337, -4.0980762113533257818, 1. ], [ 2.3660254037844303809, -5.0980762113533248936, 1. ], 
 [ 3.2320508075688669791, -5.5980762113533284463, 1. ] ], [ [ 4.4641016151377606036, -9.2148511043887992855e-15, 0. ], [ 3.4641016151377614918, -8.2156503822261583991e-15, 0. ], 
 [ 3.46410161513776238, 0.99999999999999245048, 0. ], [ 2.5980762113533248936, 1.4999999999999942268, 0. ], [ 1.7320508075688867411, 1.9999999999999955591, 0. ], 
 [ 0.8660254037844485886, 1.4999999999999964473, 0. ], [ 1.3660254037844468122, 0.63397459621555629639, 0. ], [ 0.86602540378444547997, -0.23205080756888196714, 0. ], 
 [ 1.7320508075688827443, -0.73205080756888407656, 0. ], [ 1.7320508075688803018, -1.7320508075688847427, 0. ], [ 2.7320508075688794136, -1.7320508075688874072, 0. ], 
 [ 3.23205080756888119, -0.86602540378444947677, 0. ], [ 4.0980762113533177882, -1.3660254037844521413, 0. ], [ 4.9641016151377570509, -0.86602540378445402869, 0. ], 
 [ 4.2141016151377606036, 0.24999999999999111822, 1. ], [ 3.2141016151377619359, 0.74999999999999245048, 1. ], [ 3.3725952641916530084, 1.3415063509461027103, 1. ], 
 [ 1.8235571585149963347, 1.6584936490538852993, 1. ], [ 1.3905444566227775915, 2.0915063509461058189, 1. ], [ 1.4575317547305566279, 0.97548094716166633411, 1. ], 
 [ 1.0245190528383369966, 0.54246824526944736888, 1. ], [ 1.6405444566227733727, -0.39054445662277381679, 1. ], [ 1.4820508075688823002, -0.98205080756888396554, 1. ], 
 [ 2.4820508075688798577, -1.4820508075688865191, 1. ], [ 3.0735571585149887852, -1.6405444566227782577, 1. ], [ 4.0065698604072084166, -1.0245190528383418815, 1. ], 
 [ 4.4395825622994271598, -1.4575317547305624011, 1. ], [ 4.3725952641916494557, -0.34150635094612324938, 1. ], [ 4.4641016151377606036, -9.2148511043887992855e-15, 1. ], 
 [ 3.4641016151377614918, -8.2156503822261583991e-15, 1. ], [ 3.46410161513776238, 0.99999999999999245048, 1. ], [ 2.5980762113533248936, 1.4999999999999942268, 1. ], 
 [ 1.7320508075688867411, 1.9999999999999955591, 1. ], [ 0.8660254037844485886, 1.4999999999999964473, 1. ], [ 1.3660254037844468122, 0.63397459621555629639, 1. ], 
 [ 0.86602540378444547997, -0.23205080756888196714, 1. ], [ 1.7320508075688827443, -0.73205080756888407656, 1. ], [ 1.7320508075688803018, -1.7320508075688847427, 1. ], 
 [ 2.7320508075688794136, -1.7320508075688874072, 1. ], [ 3.23205080756888119, -0.86602540378444947677, 1. ], [ 4.0980762113533177882, -1.3660254037844521413, 1. ], 
 [ 4.9641016151377570509, -0.86602540378445402869, 1. ] ], [ [ 6.1961524227066417936, 1.73205080756887142, 0. ], [ 5.1961524227066426818, 1.7320508075688709759, 0. ], 
 [ 5.1961524227066426818, 2.7320508075688718641, 0. ], [ 4.3301270189222051954, 3.2320508075688718641, 0. ], [ 3.464101615137767709, 3.7320508075688718641, 0. ], 
 [ 2.5980762113533293345, 3.2320508075688718641, 0. ], [ 3.0980762113533275581, 2.3660254037844321573, 0. ], [ 2.59807621135332667, 1.4999999999999933387, 0. ], 
 [ 3.4641016151377637122, 0.99999999999999233946, 0. ], [ 3.4641016151377614918, -8.2156503822261583991e-15, 0. ], [ 4.4641016151377606036, -9.2148511043887992855e-15, 0. ], 
 [ 4.96410161513776238, 0.86602540378442949276, 0. ], [ 5.8301270189221989781, 0.36602540378442827151, 0. ], [ 6.6961524227066382409, 0.86602540378442760538, 0. ], 
 [ 5.9461524227066417936, 1.9820508075688718641, 1. ], [ 4.9461524227066426818, 2.4820508075688718641, 1. ], [ 5.1046460717605333102, 3.0735571585149816798, 1. ], 
 [ 3.5556079660838766365, 3.3905444566227624925, 1. ], [ 3.1225952641916578933, 3.8235571585149821239, 1. ], [ 3.1895825622994378179, 2.7075317547305424171, 1. ], 
 [ 2.7565698604072181865, 2.2745190528383227857, 1. ], [ 3.3725952641916543406, 1.3415063509461022662, 1. ], [ 3.2141016151377632681, 0.74999999999999222844, 1. ], 
 [ 4.2141016151377606036, 0.2499999999999910627, 1. ], [ 4.8056079660838699752, 0.091506350946100434296, 1. ], [ 5.7386206679760896066, 0.70753175473053830924, 1. ], 
 [ 6.1716333698683083497, 0.2745190528383182893, 1. ], [ 6.1046460717605306456, 1.3905444566227576075, 1. ], [ 6.1961524227066417936, 1.73205080756887142, 1. ], 
 [ 5.1961524227066426818, 1.7320508075688709759, 1. ], [ 5.1961524227066426818, 2.7320508075688718641, 1. ], [ 4.3301270189222051954, 3.2320508075688718641, 1. ], 
 [ 3.464101615137767709, 3.7320508075688718641, 1. ], [ 2.5980762113533293345, 3.2320508075688718641, 1. ], [ 3.0980762113533275581, 2.3660254037844321573, 1. ], 
 [ 2.59807621135332667, 1.4999999999999933387, 1. ], [ 3.4641016151377637122, 0.99999999999999233946, 1. ], [ 3.4641016151377614918, -8.2156503822261583991e-15, 1. ], 
 [ 4.4641016151377606036, -9.2148511043887992855e-15, 1. ], [ 4.96410161513776238, 0.86602540378442949276, 1. ], [ 5.8301270189221989781, 0.36602540378442827151, 1. ], 
 [ 6.6961524227066382409, 0.86602540378442760538, 1. ] ], [ [ 2.4641016151377668209, 4.7320508075688767491, 0. ], [ 1.964101615137767709, 3.8660254037844361541, 0. ], 
 [ 1.0980762113533297786, 4.36602540378443571, 0. ], [ 0.23205080756889184812, 3.8660254037844348218, 0. ], [ -0.63397459621554652642, 3.3660254037844339337, 0. ], 
 [ -0.63397459621554652642, 2.3660254037844321573, 0. ], [ 0.3660254037844525854, 2.3660254037844330455, 0. ], [ 0.86602540378445169722, 1.4999999999999933387, 0. ], 
 [ 1.7320508075688898497, 1.9999999999999935607, 0. ], [ 2.598076211353327114, 1.4999999999999928946, 0. ], [ 3.0980762113533275581, 2.3660254037844321573, 0. ], 
 [ 2.5980762113533293345, 3.2320508075688718641, 0. ], [ 3.464101615137767709, 3.7320508075688718641, 0. ], [ 3.4641016151377690413, 4.7320508075688731964, 0. ], 
 [ 2.1225952641916574493, 4.6405444566227664893, 1. ], [ 1.1895825622994391502, 4.0245190528383263384, 1. ], [ 0.75656986040722040698, 4.4575317547305459698, 1. ], 
 [ -0.29246824526943671074, 3.2745190528383236739, 1. ], [ -0.88397459621554652642, 3.1160254037844330455, 1. ], [ 0.11602540378445302949, 2.6160254037844334896, 1. ], 
 [ 0.27451905283834276972, 2.0245190528383227857, 1. ], [ 1.3905444566227802561, 2.0915063509461035984, 1. ], [ 1.8235571585149994434, 1.6584936490538839671, 1. ], 
 [ 2.7565698604072181865, 2.2745190528383218975, 1. ], [ 3.1895825622994373738, 2.707531754730541973, 1. ], [ 3.1225952641916583374, 3.8235571585149816798, 1. ], 
 [ 3.714101615137767709, 3.9820508075688723082, 1. ], [ 2.7141016151377699295, 4.4820508075688731964, 1. ], [ 2.4641016151377668209, 4.7320508075688767491, 1. ], 
 [ 1.964101615137767709, 3.8660254037844361541, 1. ], [ 1.0980762113533297786, 4.36602540378443571, 1. ], [ 0.23205080756889184812, 3.8660254037844348218, 1. ], 
 [ -0.63397459621554652642, 3.3660254037844339337, 1. ], [ -0.63397459621554652642, 2.3660254037844321573, 1. ], [ 0.3660254037844525854, 2.3660254037844330455, 1. ], 
 [ 0.86602540378445169722, 1.4999999999999933387, 1. ], [ 1.7320508075688898497, 1.9999999999999935607, 1. ], [ 2.598076211353327114, 1.4999999999999928946, 1. ], 
 [ 3.0980762113533275581, 2.3660254037844321573, 1. ], [ 2.5980762113533293345, 3.2320508075688718641, 1. ], [ 3.464101615137767709, 3.7320508075688718641, 1. ], 
 [ 3.4641016151377690413, 4.7320508075688731964, 1. ] ], [ [ -1.9999999999999882316, 2.9999999999999920064, 0. ], [ -1.4999999999999873435, 2.1339745962155531878, 0. ], 
 [ -2.3660254037844250519, 1.6339745962155500791, 0. ], [ -2.3660254037844228314, 0.63397459621554785869, 0. ], [ -2.3660254037844219432, -0.36602540378445391767, 0. ], 
 [ -1.4999999999999833467, -0.86602540378445347358, 0. ], [ -0.99999999999998490097, -1.1990408665951690637e-14, 0. ], [ 1.4765966227514581988e-14, -1.0214051826551440172e-14, 0. ], 
 [ 1.4321877017664519371e-14, 0.99999999999999156231, 0. ], [ 0.86602540378445214131, 1.4999999999999935607, 0. ], [ 0.36602540378445269642, 2.3660254037844330455, 0. ], 
 [ -0.63397459621554652642, 2.3660254037844321573, 0. ], [ -0.63397459621554652642, 3.3660254037844339337, 0. ], [ -1.4999999999999844569, 3.8660254037844339337, 0. ], 
 [ -2.0915063509460969371, 2.6584936490538817466, 1. ], [ -2.0245190528383147921, 1.5424682452694407075, 1. ], [ -2.6160254037844241637, 1.3839745962155487469, 1. ], 
 [ -2.1160254037844223873, -0.11602540378445302949, 1. ], [ -2.2745190528383121276, -0.70753175473056373335, 1. ], [ -1.3415063509460942726, -0.091506350946122694268, 1. ], 
 [ -0.74999999999998467892, -0.25000000000001199041, 1. ], [ -0.24999999999998545608, 0.74999999999999067413, 1. ], [ 0.34150635094612380449, 0.90849364905388219071, 1. ], 
 [ 0.27451905283834332483, 2.0245190528383227857, 1. ], [ 0.1160254037844525854, 2.6160254037844330455, 1. ], [ -0.88397459621554619336, 3.1160254037844330455, 1. ], 
 [ -0.7254809471616563421, 3.7075317547305441934, 1. ], [ -1.6584936490538746412, 3.0915063509461040425, 1. ], [ -1.9999999999999882316, 2.9999999999999920064, 1. ], 
 [ -1.4999999999999873435, 2.1339745962155531878, 1. ], [ -2.3660254037844250519, 1.6339745962155500791, 1. ], [ -2.3660254037844228314, 0.63397459621554785869, 1. ], 
 [ -2.3660254037844219432, -0.36602540378445391767, 1. ], [ -1.4999999999999833467, -0.86602540378445347358, 1. ], [ -0.99999999999998490097, -1.1990408665951690637e-14, 1. ], 
 [ 1.4765966227514581988e-14, -1.0214051826551440172e-14, 1. ], [ 1.4321877017664519371e-14, 0.99999999999999156231, 1. ], [ 0.86602540378445214131, 1.4999999999999935607, 1. ], 
 [ 0.36602540378445269642, 2.3660254037844330455, 1. ], [ -0.63397459621554652642, 2.3660254037844321573, 1. ], [ -0.63397459621554652642, 3.3660254037844339337, 1. ], 
 [ -1.4999999999999844569, 3.8660254037844339337, 1. ] ], [ [ -4.3660254037844321573, 3.633974596215549635, 0. ], [ -3.8660254037844294928, 2.7679491924311112605, 0. ], 
 [ -4.7320508075688678673, 2.2679491924311081519, 0. ], [ -4.7320508075688652028, 1.2679491924311068196, 0. ], [ -4.7320508075688625382, 0.26794919243110548734, 0. ], 
 [ -3.8660254037844232755, -0.23205080756889318039, 0. ], [ -3.3660254037844241637, 0.63397459621554830278, 0. ], [ -2.3660254037844241637, 0.63397459621555007914, 0. ], 
 [ -2.366025403784425496, 1.6339745962155514114, 0. ], [ -1.4999999999999873435, 2.1339745962155531878, 0. ], [ -1.9999999999999882316, 2.9999999999999920064, 0. ], 
 [ -2.9999999999999884537, 2.9999999999999906741, 0. ], [ -2.9999999999999893419, 3.9999999999999920064, 0. ], [ -3.8660254037844286046, 4.4999999999999911182, 0. ], 
 [ -4.4575317547305406407, 3.2924682452694389312, 1. ], [ -4.3905444566227576075, 2.1764428414849992244, 1. ], [ -4.9820508075688669791, 2.0179491924311072637, 1. ], 
 [ -4.4820508075688634264, 0.51794919243110637552, 1. ], [ -4.6405444566227522785, -0.073557158515004772426, 1. ], [ -3.7075317547305344235, 0.54246824526943715483, 1. ], 
 [ -3.1160254037844241637, 0.38397459621554830278, 1. ], [ -2.6160254037844250519, 1.3839745962155500791, 1. ], [ -2.0245190528383156803, 1.5424682452694424839, 1. ], 
 [ -2.0915063509460973812, 2.6584936490538817466, 1. ], [ -2.2499999999999888978, 3.2499999999999920064, 1. ], [ -3.2499999999999893419, 3.7499999999999911182, 1. ], 
 [ -3.0915063509460996016, 4.341506350946101378, 1. ], [ -4.0245190528383183448, 3.7254809471616612271, 1. ], [ -4.3660254037844321573, 3.633974596215549635, 1. ], 
 [ -3.8660254037844294928, 2.7679491924311112605, 1. ], [ -4.7320508075688678673, 2.2679491924311081519, 1. ], [ -4.7320508075688652028, 1.2679491924311068196, 1. ], 
 [ -4.7320508075688625382, 0.26794919243110548734, 1. ], [ -3.8660254037844232755, -0.23205080756889318039, 1. ], [ -3.3660254037844241637, 0.63397459621554830278, 1. ], 
 [ -2.3660254037844241637, 0.63397459621555007914, 1. ], [ -2.366025403784425496, 1.6339745962155514114, 1. ], [ -1.4999999999999873435, 2.1339745962155531878, 1. ], 
 [ -1.9999999999999882316, 2.9999999999999920064, 1. ], [ -2.9999999999999884537, 2.9999999999999906741, 1. ], [ -2.9999999999999893419, 3.9999999999999920064, 1. ], 
 [ -3.8660254037844286046, 4.4999999999999911182, 1. ] ], [ [ -5.0980762113533009128, -1.0980762113533411028, 0. ], [ -4.0980762113533009128, -1.098076211353336884, 0. ], 
 [ -4.0980762113532982482, -2.0980762113533391044, 0. ], [ -3.2320508075688580973, -2.5980762113533364399, 0. ], [ -2.3660254037844179464, -3.0980762113533346636, 0. ], 
 [ -1.499999999999980016, -2.5980762113533302227, 0. ], [ -1.9999999999999822364, -1.7320508075688916261, 0. ], [ -1.4999999999999831246, -0.86602540378444925473, 0. ], 
 [ -2.3660254037844228314, -0.36602540378445103109, 0. ], [ -2.3660254037844241637, 0.63397459621555096732, 0. ], [ -3.3660254037844241637, 0.63397459621554830278, 0. ], 
 [ -3.8660254037844232755, -0.2320508075688932359, 0. ], [ -4.7320508075688625382, 0.26794919243110548734, 0. ], [ -5.5980762113533009128, -0.23205080756889762128, 0. ], 
 [ -4.8480762113533000246, -1.3480762113533404367, 1. ], [ -3.8480762113532982482, -1.8480762113533377722, 1. ], [ -4.0065698604071879885, -2.4395825622994493642, 1. ], 
 [ -2.4575317547305282062, -2.7565698604072244038, 1. ], [ -2.0245190528383076867, -3.189582562299443147, 1. ], [ -2.0915063509460911639, -2.073557158515002552, 1. ], 
 [ -1.6584936490538724208, -1.6405444566227804781, 1. ], [ -2.2745190528383125717, -0.70753175473056129086, 1. ], [ -2.1160254037844237196, -0.11602540378444953229, 1. ], 
 [ -3.1160254037844232755, 0.38397459621554863585, 1. ], [ -3.7075317547305335353, 0.54246824526943737688, 1. ], [ -4.6405444566227522785, -0.073557158515004716914, 1. ], 
 [ -5.073557158514972798, 0.35945554337721430382, 1. ], [ -5.006569860407190653, -0.7565698604072259581, 1. ], [ -5.0980762113533009128, -1.0980762113533411028, 1. ], 
 [ -4.0980762113533009128, -1.098076211353336884, 1. ], [ -4.0980762113532982482, -2.0980762113533391044, 1. ], [ -3.2320508075688580973, -2.5980762113533364399, 1. ], 
 [ -2.3660254037844179464, -3.0980762113533346636, 1. ], [ -1.499999999999980016, -2.5980762113533302227, 1. ], [ -1.9999999999999822364, -1.7320508075688916261, 1. ], 
 [ -1.4999999999999831246, -0.86602540378444925473, 1. ], [ -2.3660254037844228314, -0.36602540378445103109, 1. ], [ -2.3660254037844241637, 0.63397459621555096732, 1. ], 
 [ -3.3660254037844241637, 0.63397459621554830278, 1. ], [ -3.8660254037844232755, -0.2320508075688932359, 1. ], [ -4.7320508075688625382, 0.26794919243110548734, 1. ], 
 [ -5.5980762113533009128, -0.23205080756889762128, 1. ] ] ]


verticesOfFaces=[ [ [ 1, 2, 12 ], [ 1, 2, 15 ], [ 1, 12, 13 ], [ 1, 13, 14 ], [ 1, 29, 15 ], [ 1, 29, 28 ], [ 2, 3, 4 ], [ 2, 3, 16 ], [ 2, 9, 12 ], [ 2, 30, 15 ], 
[ 2, 30, 16 ], [ 3, 4, 17 ], [ 3, 31, 16 ], [ 3, 31, 17 ], [ 4, 5, 7 ], [ 4, 5, 18 ], [ 4, 7, 2 ], [ 4, 32, 17 ], [ 4, 32, 18 ], [ 5, 6, 7 ], [ 5, 6, 19 ], [ 5, 33, 18 ], [ 5, 33, 19 ], 
[ 6, 7, 20 ], [ 6, 34, 19 ], [ 6, 34, 20 ], [ 7, 8, 9 ], [ 7, 8, 21 ], [ 7, 35, 20 ], [ 7, 35, 21 ], [ 8, 9, 22 ], [ 8, 36, 21 ], [ 8, 36, 22 ], 
[ 9, 7, 2 ], [ 9, 10, 11 ], [ 9, 10, 23 ], [ 9, 11, 12 ], [ 9, 37, 22 ], [ 9, 37, 23 ], [ 10, 11, 24 ], [ 10, 38, 23 ], [ 10, 38, 24 ], [ 11, 12, 25 ], [ 11, 39, 24 ], [ 11, 39, 25 ], 
[ 12, 13, 26 ], [ 12, 40, 25 ], [ 12, 40, 26 ], [ 13, 14, 27 ], [ 13, 41, 26 ], [ 13, 41, 27 ], [ 14, 1, 28 ], [ 14, 42, 27 ], [ 14, 42, 28 ], 
[ 16, 17, 31 ], [ 16, 17, 32 ], [ 16, 18, 20 ], [ 16, 18, 32 ], [ 16, 20, 35 ], [ 16, 22, 35 ], [ 16, 24, 38 ], [ 16, 37, 22 ], [ 16, 38, 37 ], [ 18, 19, 33 ], [ 18, 19, 34 ], [ 18, 34, 20 ], 
[ 21, 35, 36 ], [ 24, 25, 40 ], [ 24, 30, 16 ], [ 24, 39, 25 ], [ 24, 40, 30 ], [ 26, 27, 41 ], [ 27, 28, 26 ], [ 27, 28, 42 ], [ 35, 36, 22 ], [ 37, 38, 23 ], 
[ 40, 26, 28 ], [ 40, 28, 29 ], [ 40, 29, 15 ], [ 40, 30, 15 ] ], [ [ 43, 44, 54 ], [ 43, 44, 57 ], [ 43, 54, 55 ], [ 43, 55, 56 ], [ 43, 71, 57 ], [ 43, 71, 70 ], [ 44, 45, 46 ], 
[ 44, 45, 58 ], [ 44, 51, 54 ], [ 44, 72, 57 ], [ 44, 72, 58 ], [ 45, 46, 59 ], [ 45, 73, 58 ], [ 45, 73, 59 ], [ 46, 47, 49 ], [ 46, 47, 60 ], [ 46, 49, 44 ], 
[ 46, 74, 59 ], [ 46, 74, 60 ], [ 47, 48, 49 ], [ 47, 48, 61 ], [ 47, 75, 60 ], [ 47, 75, 61 ], [ 48, 49, 62 ], [ 48, 76, 61 ], [ 48, 76, 62 ], [ 49, 50, 51 ], [ 49, 50, 63 ], [ 49, 77, 62 ], 
[ 49, 77, 63 ], [ 50, 51, 64 ], [ 50, 78, 63 ], [ 50, 78, 64 ], [ 51, 49, 44 ], [ 51, 52, 53 ], [ 51, 52, 65 ], [ 51, 53, 54 ], [ 51, 79, 64 ], [ 51, 79, 65 ], 
[ 52, 53, 66 ], [ 52, 80, 65 ], [ 52, 80, 66 ], [ 53, 54, 67 ], [ 53, 81, 66 ], [ 53, 81, 67 ], [ 54, 55, 68 ], [ 54, 82, 67 ], [ 54, 82, 68 ], [ 55, 56, 69 ], [ 55, 83, 68 ], [ 55, 83, 69 ], 
[ 56, 43, 70 ], [ 56, 84, 69 ], [ 56, 84, 70 ], [ 58, 59, 73 ], [ 58, 59, 74 ], [ 58, 60, 62 ], [ 58, 60, 74 ], [ 58, 62, 77 ], [ 58, 64, 77 ], [ 58, 66, 80 ], 
[ 58, 79, 64 ], [ 58, 80, 79 ], [ 60, 61, 75 ], [ 60, 61, 76 ], [ 60, 76, 62 ], [ 63, 77, 78 ], [ 66, 67, 82 ], [ 66, 72, 58 ], [ 66, 81, 67 ], [ 66, 82, 72 ], [ 68, 69, 83 ], [ 69, 70, 68 ], 
[ 69, 70, 84 ], [ 77, 78, 64 ], [ 79, 80, 65 ], [ 82, 68, 70 ], [ 82, 70, 71 ], [ 82, 71, 57 ], [ 82, 72, 57 ] ], 
[ [ 85, 86, 96 ], [ 85, 86, 99 ], [ 85, 96, 97 ], [ 85, 97, 98 ], [ 85, 113, 99 ], [ 85, 113, 112 ], [ 86, 87, 88 ], 
[ 86, 87, 100 ], [ 86, 93, 96 ], [ 86, 114, 99 ], [ 86, 114, 100 ], [ 87, 88, 101 ], [ 87, 115, 100 ], [ 87, 115, 101 ], [ 88, 89, 91 ], [ 88, 89, 102 ], [ 88, 91, 86 ], [ 88, 116, 101 ], 
[ 88, 116, 102 ], [ 89, 90, 91 ], [ 89, 90, 103 ], [ 89, 117, 102 ], [ 89, 117, 103 ], [ 90, 91, 104 ], [ 90, 118, 103 ], [ 90, 118, 104 ], [ 91, 92, 93 ], [ 91, 92, 105 ], 
[ 91, 119, 104 ], [ 91, 119, 105 ], [ 92, 93, 106 ], [ 92, 120, 105 ], [ 92, 120, 106 ], [ 93, 91, 86 ], [ 93, 94, 95 ], [ 93, 94, 107 ], [ 93, 95, 96 ], [ 93, 121, 106 ], [ 93, 121, 107 ], 
[ 94, 95, 108 ], [ 94, 122, 107 ], [ 94, 122, 108 ], [ 95, 96, 109 ], [ 95, 123, 108 ], [ 95, 123, 109 ], [ 96, 97, 110 ], [ 96, 124, 109 ], [ 96, 124, 110 ], [ 97, 98, 111 ], 
[ 97, 125, 110 ], [ 97, 125, 111 ], [ 98, 85, 112 ], [ 98, 126, 111 ], [ 98, 126, 112 ], [ 100, 101, 115 ], [ 100, 101, 116 ], [ 100, 102, 104 ], [ 100, 102, 116 ], [ 100, 104, 119 ], 
[ 100, 106, 119 ], [ 100, 108, 122 ], [ 100, 121, 106 ], [ 100, 122, 121 ], [ 102, 103, 117 ], [ 102, 103, 118 ], [ 102, 118, 104 ], [ 105, 119, 120 ], [ 108, 109, 124 ], 
[ 108, 114, 100 ], [ 108, 123, 109 ], [ 108, 124, 114 ], [ 110, 111, 125 ], [ 111, 112, 110 ], [ 111, 112, 126 ], [ 119, 120, 106 ], [ 121, 122, 107 ], [ 124, 110, 112 ], [ 124, 112, 113 ], 
[ 124, 113, 99 ], [ 124, 114, 99 ] ], [ [ 127, 128, 138 ], [ 127, 128, 141 ], [ 127, 138, 139 ], [ 127, 139, 140 ], [ 127, 155, 141 ], [ 127, 155, 154 ], [ 128, 129, 130 ], 
[ 128, 129, 142 ], [ 128, 135, 138 ], [ 128, 156, 141 ], [ 128, 156, 142 ], [ 129, 130, 143 ], [ 129, 157, 142 ], [ 129, 157, 143 ], [ 130, 131, 133 ], [ 130, 131, 144 ], [ 130, 133, 128 ], 
[ 130, 158, 143 ], [ 130, 158, 144 ], [ 131, 132, 133 ], [ 131, 132, 145 ], [ 131, 159, 144 ], [ 131, 159, 145 ], [ 132, 133, 146 ], [ 132, 160, 145 ], [ 132, 160, 146 ], [ 133, 134, 135 ], 
[ 133, 134, 147 ], [ 133, 161, 146 ], [ 133, 161, 147 ], [ 134, 135, 148 ], [ 134, 162, 147 ], [ 134, 162, 148 ], [ 135, 133, 128 ], [ 135, 136, 137 ], [ 135, 136, 149 ], [ 135, 137, 138 ], 
[ 135, 163, 148 ], [ 135, 163, 149 ], [ 136, 137, 150 ], [ 136, 164, 149 ], [ 136, 164, 150 ], [ 137, 138, 151 ], [ 137, 165, 150 ], [ 137, 165, 151 ], [ 138, 139, 152 ], [ 138, 166, 151 ], 
[ 138, 166, 152 ], [ 139, 140, 153 ], [ 139, 167, 152 ], [ 139, 167, 153 ], [ 140, 127, 154 ], [ 140, 168, 153 ], [ 140, 168, 154 ], [ 142, 143, 157 ], [ 142, 143, 158 ], [ 142, 144, 146 ], 
[ 142, 144, 158 ], [ 142, 146, 161 ], [ 142, 148, 161 ], [ 142, 150, 164 ], [ 142, 163, 148 ], [ 142, 164, 163 ], [ 144, 145, 159 ], [ 144, 145, 160 ], [ 144, 160, 146 ], [ 147, 161, 162 ], 
[ 150, 151, 166 ], [ 150, 156, 142 ], [ 150, 165, 151 ], [ 150, 166, 156 ], [ 152, 153, 167 ], [ 153, 154, 152 ], [ 153, 154, 168 ], [ 161, 162, 148 ], [ 163, 164, 149 ], [ 166, 152, 154 ], 
[ 166, 154, 155 ], [ 166, 155, 141 ], [ 166, 156, 141 ] ], [ [ 169, 170, 180 ], [ 169, 170, 183 ], [ 169, 180, 181 ], [ 169, 181, 182 ], [ 169, 197, 183 ], [ 169, 197, 196 ], [ 170, 171, 172 ], 
[ 170, 171, 184 ], [ 170, 177, 180 ], [ 170, 198, 183 ], [ 170, 198, 184 ], [ 171, 172, 185 ], [ 171, 199, 184 ], [ 171, 199, 185 ], [ 172, 173, 175 ], [ 172, 173, 186 ], [ 172, 175, 170 ], 
[ 172, 200, 185 ], [ 172, 200, 186 ], [ 173, 174, 175 ], [ 173, 174, 187 ], [ 173, 201, 186 ], [ 173, 201, 187 ], [ 174, 175, 188 ], [ 174, 202, 187 ], [ 174, 202, 188 ], [ 175, 176, 177 ], 
[ 175, 176, 189 ], [ 175, 203, 188 ], [ 175, 203, 189 ], [ 176, 177, 190 ], [ 176, 204, 189 ], [ 176, 204, 190 ], [ 177, 175, 170 ], [ 177, 178, 179 ], [ 177, 178, 191 ], [ 177, 179, 180 ], 
[ 177, 205, 190 ], [ 177, 205, 191 ], [ 178, 179, 192 ], [ 178, 206, 191 ], [ 178, 206, 192 ], [ 179, 180, 193 ], [ 179, 207, 192 ], [ 179, 207, 193 ], [ 180, 181, 194 ], [ 180, 208, 193 ], 
[ 180, 208, 194 ], [ 181, 182, 195 ], [ 181, 209, 194 ], [ 181, 209, 195 ], [ 182, 169, 196 ], [ 182, 210, 195 ], [ 182, 210, 196 ], [ 184, 185, 199 ], [ 184, 185, 200 ], [ 184, 186, 188 ], 
[ 184, 186, 200 ], [ 184, 188, 203 ], [ 184, 190, 203 ], [ 184, 192, 206 ], [ 184, 205, 190 ], [ 184, 206, 205 ], [ 186, 187, 201 ], [ 186, 187, 202 ], [ 186, 202, 188 ], [ 189, 203, 204 ], 
[ 192, 193, 208 ], [ 192, 198, 184 ], [ 192, 207, 193 ], [ 192, 208, 198 ], [ 194, 195, 209 ], [ 195, 196, 194 ], [ 195, 196, 210 ], [ 203, 204, 190 ], [ 205, 206, 191 ], [ 208, 194, 196 ], 
[ 208, 196, 197 ], [ 208, 197, 183 ], [ 208, 198, 183 ] ], [ [ 211, 212, 222 ], [ 211, 212, 225 ], [ 211, 222, 223 ], [ 211, 223, 224 ], [ 211, 239, 225 ], [ 211, 239, 238 ], [ 212, 213, 214 ], 
[ 212, 213, 226 ], [ 212, 219, 222 ], [ 212, 240, 225 ], [ 212, 240, 226 ], [ 213, 214, 227 ], [ 213, 241, 226 ], [ 213, 241, 227 ], [ 214, 215, 217 ], [ 214, 215, 228 ], [ 214, 217, 212 ], 
[ 214, 242, 227 ], [ 214, 242, 228 ], [ 215, 216, 217 ], [ 215, 216, 229 ], [ 215, 243, 228 ], [ 215, 243, 229 ], [ 216, 217, 230 ], [ 216, 244, 229 ], [ 216, 244, 230 ], [ 217, 218, 219 ], 
[ 217, 218, 231 ], [ 217, 245, 230 ], [ 217, 245, 231 ], [ 218, 219, 232 ], [ 218, 246, 231 ], [ 218, 246, 232 ], [ 219, 217, 212 ], [ 219, 220, 221 ], [ 219, 220, 233 ], [ 219, 221, 222 ], 
[ 219, 247, 232 ], [ 219, 247, 233 ], [ 220, 221, 234 ], [ 220, 248, 233 ], [ 220, 248, 234 ], [ 221, 222, 235 ], [ 221, 249, 234 ], [ 221, 249, 235 ], [ 222, 223, 236 ], [ 222, 250, 235 ], 
[ 222, 250, 236 ], [ 223, 224, 237 ], [ 223, 251, 236 ], [ 223, 251, 237 ], [ 224, 211, 238 ], [ 224, 252, 237 ], [ 224, 252, 238 ], [ 226, 227, 241 ], [ 226, 227, 242 ], [ 226, 228, 230 ], 
[ 226, 228, 242 ], [ 226, 230, 245 ], [ 226, 232, 245 ], [ 226, 234, 248 ], [ 226, 247, 232 ], [ 226, 248, 247 ], [ 228, 229, 243 ], [ 228, 229, 244 ], [ 228, 244, 230 ], [ 231, 245, 246 ], 
[ 234, 235, 250 ], [ 234, 240, 226 ], [ 234, 249, 235 ], [ 234, 250, 240 ], [ 236, 237, 251 ], [ 237, 238, 236 ], [ 237, 238, 252 ], [ 245, 246, 232 ], [ 247, 248, 233 ], [ 250, 236, 238 ], 
[ 250, 238, 239 ], [ 250, 239, 225 ], [ 250, 240, 225 ] ], [ [ 253, 254, 264 ], [ 253, 254, 267 ], [ 253, 264, 265 ], [ 253, 265, 266 ], [ 253, 281, 267 ], [ 253, 281, 280 ], [ 254, 255, 256 ], 
[ 254, 255, 268 ], [ 254, 261, 264 ], [ 254, 282, 267 ], [ 254, 282, 268 ], [ 255, 256, 269 ], [ 255, 283, 268 ], [ 255, 283, 269 ], [ 256, 257, 259 ], [ 256, 257, 270 ], [ 256, 259, 254 ], 
[ 256, 284, 269 ], [ 256, 284, 270 ], [ 257, 258, 259 ], [ 257, 258, 271 ], [ 257, 285, 270 ], [ 257, 285, 271 ], [ 258, 259, 272 ], [ 258, 286, 271 ], [ 258, 286, 272 ], [ 259, 260, 261 ], 
[ 259, 260, 273 ], [ 259, 287, 272 ], [ 259, 287, 273 ], [ 260, 261, 274 ], [ 260, 288, 273 ], [ 260, 288, 274 ], [ 261, 259, 254 ], [ 261, 262, 263 ], [ 261, 262, 275 ], [ 261, 263, 264 ], 
[ 261, 289, 274 ], [ 261, 289, 275 ], [ 262, 263, 276 ], [ 262, 290, 275 ], [ 262, 290, 276 ], [ 263, 264, 277 ], [ 263, 291, 276 ], [ 263, 291, 277 ], [ 264, 265, 278 ], [ 264, 292, 277 ], 
[ 264, 292, 278 ], [ 265, 266, 279 ], [ 265, 293, 278 ], [ 265, 293, 279 ], [ 266, 253, 280 ], [ 266, 294, 279 ], [ 266, 294, 280 ], [ 268, 269, 283 ], [ 268, 269, 284 ], [ 268, 270, 272 ], 
[ 268, 270, 284 ], [ 268, 272, 287 ], [ 268, 274, 287 ], [ 268, 276, 290 ], [ 268, 289, 274 ], [ 268, 290, 289 ], [ 270, 271, 285 ], [ 270, 271, 286 ], [ 270, 286, 272 ], [ 273, 287, 288 ], 
[ 276, 277, 292 ], [ 276, 282, 268 ], [ 276, 291, 277 ], [ 276, 292, 282 ], [ 278, 279, 293 ], [ 279, 280, 278 ], [ 279, 280, 294 ], [ 287, 288, 274 ], [ 289, 290, 275 ], [ 292, 278, 280 ], 
[ 292, 280, 281 ], [ 292, 281, 267 ], [ 292, 282, 267 ] ], [ [ 295, 296, 306 ], [ 295, 296, 309 ], [ 295, 306, 307 ], [ 295, 307, 308 ], [ 295, 323, 309 ], [ 295, 323, 322 ], [ 296, 297, 298 ], 
[ 296, 297, 310 ], [ 296, 303, 306 ], [ 296, 324, 309 ], [ 296, 324, 310 ], [ 297, 298, 311 ], [ 297, 325, 310 ], [ 297, 325, 311 ], [ 298, 299, 301 ], [ 298, 299, 312 ], [ 298, 301, 296 ], 
[ 298, 326, 311 ], [ 298, 326, 312 ], [ 299, 300, 301 ], [ 299, 300, 313 ], [ 299, 327, 312 ], [ 299, 327, 313 ], [ 300, 301, 314 ], [ 300, 328, 313 ], [ 300, 328, 314 ], [ 301, 302, 303 ], 
[ 301, 302, 315 ], [ 301, 329, 314 ], [ 301, 329, 315 ], [ 302, 303, 316 ], [ 302, 330, 315 ], [ 302, 330, 316 ], [ 303, 301, 296 ], [ 303, 304, 305 ], [ 303, 304, 317 ], [ 303, 305, 306 ], 
[ 303, 331, 316 ], [ 303, 331, 317 ], [ 304, 305, 318 ], [ 304, 332, 317 ], [ 304, 332, 318 ], [ 305, 306, 319 ], [ 305, 333, 318 ], [ 305, 333, 319 ], [ 306, 307, 320 ], [ 306, 334, 319 ], 
[ 306, 334, 320 ], [ 307, 308, 321 ], [ 307, 335, 320 ], [ 307, 335, 321 ], [ 308, 295, 322 ], [ 308, 336, 321 ], [ 308, 336, 322 ], [ 310, 311, 325 ], [ 310, 311, 326 ], [ 310, 312, 314 ], 
[ 310, 312, 326 ], [ 310, 314, 329 ], [ 310, 316, 329 ], [ 310, 318, 332 ], [ 310, 331, 316 ], [ 310, 332, 331 ], [ 312, 313, 327 ], [ 312, 313, 328 ], [ 312, 328, 314 ], [ 315, 329, 330 ], 
[ 318, 319, 334 ], [ 318, 324, 310 ], [ 318, 333, 319 ], [ 318, 334, 324 ], [ 320, 321, 335 ], [ 321, 322, 320 ], [ 321, 322, 336 ], [ 329, 330, 316 ], [ 331, 332, 317 ], [ 334, 320, 322 ], 
[ 334, 322, 323 ], [ 334, 323, 309 ], [ 334, 324, 309 ] ], [ [ 337, 338, 348 ], [ 337, 338, 351 ], [ 337, 348, 349 ], [ 337, 349, 350 ], [ 337, 365, 351 ], [ 337, 365, 364 ], [ 338, 339, 340 ], 
[ 338, 339, 352 ], [ 338, 345, 348 ], [ 338, 366, 351 ], [ 338, 366, 352 ], [ 339, 340, 353 ], [ 339, 367, 352 ], [ 339, 367, 353 ], [ 340, 341, 343 ], [ 340, 341, 354 ], [ 340, 343, 338 ], 
[ 340, 368, 353 ], [ 340, 368, 354 ], [ 341, 342, 343 ], [ 341, 342, 355 ], [ 341, 369, 354 ], [ 341, 369, 355 ], [ 342, 343, 356 ], [ 342, 370, 355 ], [ 342, 370, 356 ], [ 343, 344, 345 ], 
[ 343, 344, 357 ], [ 343, 371, 356 ], [ 343, 371, 357 ], [ 344, 345, 358 ], [ 344, 372, 357 ], [ 344, 372, 358 ], [ 345, 343, 338 ], [ 345, 346, 347 ], [ 345, 346, 359 ], [ 345, 347, 348 ], 
[ 345, 373, 358 ], [ 345, 373, 359 ], [ 346, 347, 360 ], [ 346, 374, 359 ], [ 346, 374, 360 ], [ 347, 348, 361 ], [ 347, 375, 360 ], [ 347, 375, 361 ], [ 348, 349, 362 ], [ 348, 376, 361 ], 
[ 348, 376, 362 ], [ 349, 350, 363 ], [ 349, 377, 362 ], [ 349, 377, 363 ], [ 350, 337, 364 ], [ 350, 378, 363 ], [ 350, 378, 364 ], [ 352, 353, 367 ], [ 352, 353, 368 ], [ 352, 354, 356 ], 
[ 352, 354, 368 ], [ 352, 356, 371 ], [ 352, 358, 371 ], [ 352, 360, 374 ], [ 352, 373, 358 ], [ 352, 374, 373 ], [ 354, 355, 369 ], [ 354, 355, 370 ], [ 354, 370, 356 ], [ 357, 371, 372 ], 
[ 360, 361, 376 ], [ 360, 366, 352 ], [ 360, 375, 361 ], [ 360, 376, 366 ], [ 362, 363, 377 ], [ 363, 364, 362 ], [ 363, 364, 378 ], [ 371, 372, 358 ], [ 373, 374, 359 ], [ 376, 362, 364 ], 
[ 376, 364, 365 ], [ 376, 365, 351 ], [ 376, 366, 351 ] ] ]

for block in verticesOfFaces
    x = min(vcat(block...)...)
    for f in block
        f .= f .- x .+ 1
    end
end

ass = map(i -> Polyhedron(verts = coordinates[i], facets = verticesOfFaces[i]), eachindex(verticesOfFaces))