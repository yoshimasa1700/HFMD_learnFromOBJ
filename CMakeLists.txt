# CMakeのバージョン
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#SET(CMAKE_CXX_COMPILER "clang++")

# リリースビルド用オプション
SET(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3 -pg -fpermissive")
SET(CMAKE_C_FLAGS_RELEASE "-fpermissive")

# デバッグビルド用オプション
SET(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -pg -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_MINSIZEREL "")
SET (CMAKE_CXX_FLAGS "-pg")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pg")

# ビルドのタイプはリリース
#SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_BUILD_TYPE Debug)
# 環境変数の登録


# プロジェクト名
PROJECT(objectPoseEstimation)

# 実行に必要なライブラリのパッケージを見つける
FIND_PACKAGE(OpenCV)
FIND_PACKAGE(Boost)
find_package(OpenMP)
find_package(OpenGL)
find_package(GLUT)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# バージョン番号
# SET(serial "1.0.0")

# ビルドするためにインクルードするヘッダーファイルの場所
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")
INCLUDE_DIRECTORIES("/usr/local/include")
INCLUDE_DIRECTORIES(${OPENCV_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
include_directories(${OPENGL_INCLUDE_DIRS})


# ライブラリの場所
LINK_DIRECTORIES("/usr/local/lib")

#SET(RANDOM_FOREST src/CRForest.cpp src/util.cpp src/CRTree.cpp source/HoG.cpp source/CDataset.cpp source/CConfig.cpp source/globjloader.cpp)

# 実行ファイルを生成
ADD_EXECUTABLE(objectPoseEstimation src/objectPoseEstimation.cpp)# ${RANDOM_FOREST})
ADD_EXECUTABLE(learning src/learning.cpp)# ${RANDOM_FOREST})

# ライブラリとリンクする
TARGET_LINK_LIBRARIES(objectPoseEstimation HFMD_core ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${OPENGL_LIBRARIES} )
TARGET_LINK_LIBRARIES(learning HFMD_core ${OPENGL_LIBRARIES} )#${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} glm png jpeg HFMD_core)
# バージョン番号の適用
# SET_TARGET_PROPERTIES(objectPoseEstimation PROPERTIES VERSION ${serial})
# SET_TARGET_PROPERTIES(learning PROPERTIES VERSION ${serial})
