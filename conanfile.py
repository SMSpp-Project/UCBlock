from conans import ConanFile, CMake, tools


class UcblockConan(ConanFile):
    name = "ucblock"
    version = "0.6.2"
    description = "Unit Commitment Blocks for SMS++"
    topics = ("conan", "smspp", "ucblock")
    url = "https://gitlab.com/smspp/ucblock"
    homepage = "https://gitlab.com/smspp/ucblock"
    license = "GPL-3.0-only"
    generators = "cmake"

    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    requires = (
        "smspp/0.5.2@smspp/testing"
    )

    exports_sources = [
        "CMakeLists.txt",
        "src/*",
        "include/*",
        "cmake/*",
        "tools/*"
    ]

    def source(self):
        tools.replace_in_file(
            "CMakeLists.txt",
            '''LANGUAGES C CXX)''',
            '''LANGUAGES C CXX)\n''' +
            '''include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)\n''' +
            '''conan_basic_setup()'''
        )

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.definitions["BUILD_TESTING"] = False
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy(pattern="LICENSE", dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.includedirs = ["include", "include/SMS++"]
        self.cpp_info.libs = ["UCBlock"]
