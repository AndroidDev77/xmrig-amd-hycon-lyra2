/* XMRig
 * Copyright 2010      Jeff Garzik <jgarzik@pobox.com>
 * Copyright 2012-2014 pooler      <pooler@litecoinpool.org>
 * Copyright 2014      Lucas Jones <https://github.com/lucasjones>
 * Copyright 2014-2016 Wolf9466    <https://github.com/OhGodAPet>
 * Copyright 2016      Jay D Dee   <jayddee246@gmail.com>
 * Copyright 2017-2018 XMR-Stak    <https://github.com/fireice-uk>, <https://github.com/psychocrypt>
 * Copyright 2018      Lee Clagett <https://github.com/vtnerd>
 * Copyright 2018-2019 SChernykh   <https://github.com/SChernykh>
 * Copyright 2016-2019 XMRig       <https://github.com/xmrig>, <support@xmrig.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>
#include <regex>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <inttypes.h>


#include "amd/OclCache.h"
#include "amd/OclError.h"
#include "amd/OclGPU.h"
#include "amd/OclLib.h"
#include "amd/OclLyra2_gen.h"
#include "common/log/Log.h"
#include "common/utils/timestamp.h"
#include "core/Config.h"
//#include "crypto/CryptoNight_constants.h"
#include "crypto/lyra2.h"


constexpr const char *kSetKernelArgErr = "Error %s when calling clSetKernelArg for kernel %d, argument %d.";


inline static const char *err_to_str(cl_int ret)
{
    return OclError::toString(ret);
}


inline static bool setKernelArgFromExtraBuffers(GpuContext *ctx, size_t kernel, cl_uint argument, size_t offset)
{
    cl_int ret;
    if ((ret = OclLib::setKernelArg(ctx->Kernels[kernel], argument, sizeof(cl_mem), ctx->ExtraBuffers + offset)) != CL_SUCCESS) {
        LOG_ERR(kSetKernelArgErr, err_to_str(ret), kernel, argument);
        return false;
    }

    return true;
}


inline static int cn0KernelOffset(xmrig::Variant variant)
{
#   ifndef XMRIG_NO_CN_GPU
    if (variant == xmrig::VARIANT_GPU) {
        return 13;
    }
#   endif

    return 0;
}

inline static int cn1KernelOffset(xmrig::Variant variant)
{
    switch (variant) {
    
    case xmrig::VARIANT_XHV:
        return 1;

    case xmrig::VARIANT_1:
    case xmrig::VARIANT_XTL:
    case xmrig::VARIANT_RTO:
        return 7;

    case xmrig::VARIANT_MSR:
        return 8;

    case xmrig::VARIANT_XAO:
        return 9;

    case xmrig::VARIANT_TUBE:
        return 10;

    case xmrig::VARIANT_2:
    case xmrig::VARIANT_TRTL:
        return 11;

    case xmrig::VARIANT_HALF:
        return 12;

    // 13, 14 reserved for cn/gpu cn0

#   ifndef XMRIG_NO_CN_GPU
    case xmrig::VARIANT_GPU:
        return 15;
#   endif

    // 16 reserved for cn/gpu cn2

    case xmrig::VARIANT_RWZ:
        return 17;

    case xmrig::VARIANT_ZLS:
        return 18;

    case xmrig::VARIANT_DOUBLE:
        return 19;
	case xmrig::VARIANT_0:
    case xmrig::VARIANT_WOW:
    case xmrig::VARIANT_4:
        return 20;

    default:
        break;
    }

    assert(false);

    return 0;
}

inline static int cn2KernelOffset(xmrig::Variant variant)
{
#   ifndef XMRIG_NO_CN_GPU
    if (variant == xmrig::VARIANT_GPU) {
        return 16;
    }
#   endif

    return 2;
}


static void printGPU(int index, GpuContext *ctx, xmrig::Config *config)
{
    const size_t memSize             = LYRA2_MEMSIZE * ctx->rawIntensity;
    constexpr const size_t byteToGiB = 1024u * 1024u * 1024u;
    size_t maximumWorkSize           = 0;

    OclLib::getDeviceInfo(ctx->DeviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &maximumWorkSize);

    ctx->name         = OclLib::getDeviceName(ctx->DeviceID);
    ctx->board        = OclLib::getDeviceBoardName(ctx->DeviceID);
    ctx->computeUnits = OclLib::getDeviceMaxComputeUnits(ctx->DeviceID);

    OclLib::getDeviceInfo(ctx->DeviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(size_t), &ctx->freeMem);
    OclLib::getDeviceInfo(ctx->DeviceID, CL_DEVICE_GLOBAL_MEM_SIZE,    sizeof(size_t), &ctx->globalMem);

    ctx->freeMem = std::min(ctx->freeMem, ctx->globalMem);

    if (ctx->name == ctx->board) {
        LOG_INFO(config->isColors() ? WHITE_BOLD("#%02d") ", GPU " WHITE_BOLD("#%02zu") " " GREEN_BOLD("%s") ", i:" WHITE_BOLD("%zu") " " GRAY("(%zu/%zu)")
                                      ", si:" WHITE_BOLD("%d/%d") ", u:" WHITE_BOLD("%d") ", cu:" WHITE_BOLD("%d")
                                    : "#%02d, GPU #%02zu (%s), i:%zu (%zu/%zu), si:%d/%d, u:%d, cu:%d",
                 index, ctx->deviceIdx, ctx->board.data(), ctx->rawIntensity, ctx->workSize, maximumWorkSize, ctx->stridedIndex, ctx->memChunk, ctx->unrollFactor, ctx->computeUnits);
    }
    else {
        LOG_INFO(config->isColors() ? WHITE_BOLD("#%02d") ", GPU " WHITE_BOLD("#%02zu") " " GREEN_BOLD("%s") " (" CYAN_BOLD("%s") "), i:" WHITE_BOLD("%zu") " " GRAY("(%zu/%zu)")
                                      ", si:" WHITE_BOLD("%d/%d") ", u:" WHITE_BOLD("%d") ", cu:" WHITE_BOLD("%d")
                                    : "#%02d, GPU #%02zu %s (%s), i:%zu (%zu/%zu), si:%d/%d, u:%d, cu:%d",
                 index, ctx->deviceIdx, ctx->board.data(), ctx->name.data(), ctx->rawIntensity, ctx->workSize, maximumWorkSize, ctx->stridedIndex, ctx->memChunk, ctx->unrollFactor, ctx->computeUnits);
    }

    LOG_INFO(config->isColors() ? "             " CYAN("%1.2f/%1.2f/%1.0f") " GB"
                                : "             %1.2f/%1.2f/%1.0f GB",
             static_cast<double>(memSize) / byteToGiB, static_cast<double>(ctx->freeMem) / byteToGiB, static_cast<double>(ctx->globalMem) / byteToGiB);
}


size_t InitOpenCLGpu(int index, cl_context opencl_ctx, GpuContext* ctx, const char* source_code, xmrig::Config *config)
{
    ctx->opencl_ctx = opencl_ctx;

    printGPU(index, ctx, config);

    cl_int ret;
    ctx->CommandQueues = OclLib::createCommandQueue(opencl_ctx, ctx->DeviceID, &ret);
    if (ret != CL_SUCCESS) {
        return OCL_ERR_API;
    }

    ctx->InputBuffer = OclLib::createBuffer(opencl_ctx, CL_MEM_READ_ONLY, 512, nullptr, &ret);
    if (ret != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clCreateBuffer to create input buffer.", err_to_str(ret));
        return OCL_ERR_API;
    }

    /*size_t g_thd = ctx->rawIntensity;
    ctx->ExtraBuffers[0] = OclLib::createBuffer(opencl_ctx, CL_MEM_READ_WRITE, LYRA2_MEMSIZE * g_thd, nullptr, &ret);
    if (ret != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clCreateBuffer to create hash scratchpads buffer.", err_to_str(ret));
        return OCL_ERR_API;
    }
	//Lyra2 Input
    ctx->ExtraBuffers[1] = OclLib::createBuffer(opencl_ctx, CL_MEM_READ_WRITE, 200 * g_thd, nullptr, &ret);
    if(ret != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clCreateBuffer to create hash states buffer.", err_to_str(ret));
        return OCL_ERR_API;
    }*/


    // Assume we may find up to 0xFF nonces in one run - it's reasonable
    ctx->OutputBuffer = OclLib::createBuffer(opencl_ctx, CL_MEM_READ_WRITE, sizeof(cl_uint) * 2, nullptr, &ret);
    if (ret != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clCreateBuffer to create output buffer.", err_to_str(ret));
        return OCL_ERR_API;
    }

    OclCache cache(index, opencl_ctx, ctx, source_code, config);
    if (!cache.load()) {
        return OCL_ERR_API;
    }

    const char *KernelNames[] = {
        "lyra2",// "JH",// "Blake",
        //"Blake", "Groestl", "JH", "Skein",
        //"cn1_monero", "cn1_msr", "cn1_xao", "cn1_tube", "cn1_v2_monero", "cn1_v2_half",
#       //ifndef XMRIG_NO_CN_GPU
        //"cn0_cn_gpu", "cn00_cn_gpu", "cn1_cn_gpu", "cn2_cn_gpu",
#       //else
        //"", "", "", "",
#       //endif
       // "cn1_v2_rwz", "cn1_v2_zls", "cn1_v2_double","lyra2",

        nullptr
    };
    for (int i = 0; KernelNames[i]; ++i) {
        if (!KernelNames[i][0]) {
            continue;
        }

        ctx->Kernels[i] = OclLib::createKernel(ctx->Program, KernelNames[i], &ret);
        if (ret != CL_SUCCESS) {
            return OCL_ERR_API;
        }
    }

    ctx->Nonce = 0;
    return 0;
}


int OclGPU::findPlatformIdx(xmrig::Config *config)
{
    assert(config->vendor() > xmrig::OCL_VENDOR_MANUAL);

#   if !defined(__APPLE__)
    char name[256]  = { 0 };
    const int index = findPlatformIdx(config->vendor(), name, sizeof name);

    if (index >= 0) {
        LOG_INFO(config->isColors() ? GREEN_BOLD("found ") WHITE_BOLD("%s") " platform index: " WHITE_BOLD("%d") ", name: " WHITE_BOLD("%s")
                                    : "found %s platform index: %d, name: %s",
                 config->vendorName(config->vendor()), index , name);
    }

    return index;
#   else
    return 0;
#   endif
}


std::vector<GpuContext> OclGPU::getDevices(xmrig::Config *config)
{
    const size_t platformIndex            = static_cast<size_t>(config->platformIndex());
    std::vector<cl_platform_id> platforms = OclLib::getPlatformIDs();
    std::vector<GpuContext> ctxVec;

    cl_uint num_devices = 0;
    OclLib::getDeviceIDs(platforms[platformIndex], CL_DEVICE_TYPE_GPU, 0, nullptr, &num_devices);
    if (num_devices == 0) {
        return ctxVec;
    }

    cl_device_id *device_list = new cl_device_id[num_devices];
    OclLib::getDeviceIDs(platforms[platformIndex], CL_DEVICE_TYPE_GPU, num_devices, device_list, nullptr);

    for (cl_uint i = 0; i < num_devices; i++) {
        GpuContext ctx;
        ctx.deviceIdx    = i;
        ctx.platformIdx  = platformIndex;
        ctx.DeviceID     = device_list[i];
        ctx.computeUnits = OclLib::getDeviceMaxComputeUnits(ctx.DeviceID);
        ctx.vendor       = OclLib::getDeviceVendor(ctx.DeviceID);

        if (ctx.vendor == xmrig::OCL_VENDOR_UNKNOWN) {
            continue;
        }

        OclLib::getDeviceInfo(ctx.DeviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(size_t), &ctx.freeMem);
        OclLib::getDeviceInfo(ctx.DeviceID, CL_DEVICE_GLOBAL_MEM_SIZE,    sizeof(size_t), &ctx.globalMem);
        // if environment variable GPU_SINGLE_ALLOC_PERCENT is not set we can not allocate the full memory
        ctx.freeMem = std::min(ctx.freeMem, ctx.globalMem);

        ctx.board = OclLib::getDeviceBoardName(ctx.DeviceID);
        ctx.name  = OclLib::getDeviceName(ctx.DeviceID);

        if (ctx.board == ctx.name) {
            LOG_INFO(config->isColors() ? GREEN_BOLD("found") " OpenCL GPU: " GREEN_BOLD("%s") ", cu: " WHITE_BOLD("%d")
                                        : "found OpenCL GPU: %s, cu:",
                     ctx.name.data(), ctx.computeUnits);
        }
        else {
            LOG_INFO(config->isColors() ? GREEN_BOLD("found") " OpenCL GPU: " GREEN_BOLD("%s") " ("  CYAN_BOLD("%s") "), cu: " WHITE_BOLD("%d")
                                        : "found OpenCL GPU: %s (%s), cu:",
                     ctx.board.data(), ctx.name.data(), ctx.computeUnits);
        }

        ctxVec.push_back(ctx);
    }


    delete [] device_list;

    return ctxVec;
}


int OclGPU::findPlatformIdx(xmrig::OclVendor vendor, char *name, size_t nameSize)
{
#   if !defined(__APPLE__)
    std::vector<cl_platform_id> platforms = OclLib::getPlatformIDs();
    if (platforms.empty()) {
        return -1;
    }

    for (uint32_t i = 0; i < platforms.size(); i++) {
        OclLib::getPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, nameSize, name, nullptr);

        switch (vendor) {
        case xmrig::OCL_VENDOR_AMD:
            if (strstr(name, "Advanced Micro Devices") != nullptr) {
                return i;
            }
            break;

        case xmrig::OCL_VENDOR_NVIDIA:
            if (strstr(name, "NVIDIA") != nullptr) {
                return i;
            }
            break;

        case xmrig::OCL_VENDOR_INTEL:
            if (strstr(name, "Intel") != nullptr) {
                return i;
            }
            break;

        default:
            break;
        }
    }

    return -1;
#   else
    return 0;
#   endif
}


void printPlatforms()
{
    std::vector<cl_platform_id> platforms = OclLib::getPlatformIDs();

    char buf[128] = { 0 };

    for (size_t i = 0; i < platforms.size(); i++) {
        if (OclLib::getPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(buf), buf, nullptr) != CL_SUCCESS) {
            continue;
        }

        printf("#%zu: %s\n", i, buf);
    }
}


// RequestedDeviceIdxs is a list of OpenCL device indexes
// NumDevicesRequested is number of devices in RequestedDeviceIdxs list
// Returns 0 on success, -1 on stupid params, -2 on OpenCL API error
size_t InitOpenCL(const std::vector<GpuContext *> &contexts, xmrig::Config *config, cl_context *opencl_ctx)
{
    const size_t num_gpus                       = contexts.size();
    const size_t platform_idx                   = static_cast<size_t>(config->platformIndex());
    const std::vector<cl_platform_id> platforms = OclLib::getPlatformIDs();

    if (platforms.empty()) {
        return OCL_ERR_API;
    }

    if (platforms.size() <= platform_idx) {
        return OCL_ERR_BAD_PARAMS;
    }

    cl_int ret;
    cl_uint entries;
    if ((ret = OclLib::getDeviceIDs(platforms[platform_idx], CL_DEVICE_TYPE_GPU, 0, nullptr, &entries)) != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clGetDeviceIDs for number of devices.", err_to_str(ret));
        return OCL_ERR_API;
    }

    // Same as the platform index sanity check, except we must check all requested device indexes
    // TODO remove duplicated checks, see xmrig::Config::filter Threads()
    for (size_t i = 0; i < num_gpus; ++i) {
        if (entries <= contexts[i]->deviceIdx) {
            LOG_ERR("Selected OpenCL device index %lu doesn't exist.\n", contexts[i]->deviceIdx);
            return OCL_ERR_BAD_PARAMS;
        }
    }

#   ifdef __GNUC__
    cl_device_id DeviceIDList[entries];
#   else
    cl_device_id* DeviceIDList = (cl_device_id*)_alloca(entries * sizeof(cl_device_id));
#   endif

    if ((ret = OclLib::getDeviceIDs(platforms[platform_idx], CL_DEVICE_TYPE_GPU, entries, DeviceIDList, nullptr)) != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clGetDeviceIDs for device ID information.", err_to_str(ret));
        return OCL_ERR_API;
    }

    // Indexes sanity checked above
#   ifdef __GNUC__
    cl_device_id TempDeviceList[num_gpus];
#   else
    cl_device_id* TempDeviceList = (cl_device_id*)_alloca(num_gpus * sizeof(cl_device_id));
#   endif

    for (size_t i = 0; i < num_gpus; ++i) {
        TempDeviceList[i] = DeviceIDList[contexts[i]->deviceIdx];
    }

    *opencl_ctx = OclLib::createContext(nullptr, num_gpus, TempDeviceList, nullptr, nullptr, &ret);

    for (size_t i = 0; i < num_gpus; ++i) {
        contexts[i]->threadIdx   = i;
        contexts[i]->opencl_ctx  = *opencl_ctx;
        contexts[i]->platformIdx = platform_idx;
        contexts[i]->DeviceID    = DeviceIDList[contexts[i]->deviceIdx];
        OclCache::get_device_string(contexts[i]->platformIdx, contexts[i]->DeviceID, contexts[i]->DeviceString);
        contexts[i]->amdDriverMajorVersion = OclCache::amdDriverMajorVersion(contexts[0]);
    }

    if (ret != CL_SUCCESS) {
        return OCL_ERR_API;
    }

    /*const char *cryptonightCL =
            #include "./opencl/cryptonight.cl"
    ;
    const char *cryptonightCL2 =
            #include "./opencl/cryptonight2.cl"
    ;
    const char *blake256CL =
            #include "./opencl/blake256.cl"
    ;
    const char *groestl256CL =
            #include "./opencl/groestl256.cl"
    ;
    const char *jhCL =
            #include "./opencl/jh.cl"
    ;
    const char *wolfAesCL =
            #include "./opencl/wolf-aes.cl"
    ;
    const char *wolfSkeinCL =
            #include "./opencl/wolf-skein.cl"
    ;
    const char *fastIntMathV2CL =
        #include "./opencl/fast_int_math_v2.cl"
    ;
    const char *fastDivHeavyCL =
        #include "./opencl/fast_div_heavy.cl"
    ;
    const char *cryptonight_gpu =
        #include "./opencl/cryptonight_gpu.cl"
    ;*/

	const char* lyra2CL =
        #include "./opencl/lyra2.cl"
	;
	//const char* jhCL =
    //#include "./opencl/jh.cl"
	//;

    std::string source_code(lyra2CL);
    //source_code.append(cryptonightCL2);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_WOLF_AES"),         wolfAesCL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_WOLF_SKEIN"),       wolfSkeinCL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_JH"),               jhCL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_BLAKE256"),         blake256CL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_GROESTL256"),       groestl256CL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_FAST_INT_MATH_V2"), fastIntMathV2CL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_FAST_DIV_HEAVY"),   fastDivHeavyCL);
    //source_code = std::regex_replace(source_code, std::regex("XMRIG_INCLUDE_CN_GPU"),           cryptonight_gpu);
	

    for (size_t i = 0; i < num_gpus; ++i) {
        if (contexts[i]->stridedIndex == 2 && (contexts[i]->rawIntensity % contexts[i]->workSize) != 0) {
            const size_t reduced_intensity = (contexts[i]->rawIntensity / contexts[i]->workSize) * contexts[i]->workSize;
            contexts[i]->rawIntensity = reduced_intensity;

            LOG_WARN("AMD GPU #%zu: intensity is not a multiple of 'worksize', auto reduce intensity to %zu", contexts[i]->deviceIdx, reduced_intensity);
        }

        if (contexts[i]->rawIntensity % contexts[i]->workSize == 0) {
            contexts[i]->compMode = 0;
        }

        if ((ret = InitOpenCLGpu(i, *opencl_ctx, contexts[i], source_code.c_str(), config)) != OCL_ERR_SUCCESS) {
            return ret;
        }
    }

    return OCL_ERR_SUCCESS;
}
// target = password or hash uint64
// input = blob  char*
// input_len = size(input) uint64
size_t XMRSetJob(GpuContext *ctx, uint8_t *input, size_t input_len, uint64_t target, xmrig::Variant variant, uint64_t height)
{
    cl_int ret;
    
    cl_uint numThreads = ctx->rawIntensity;

    if ((ret = OclLib::enqueueWriteBuffer(ctx->CommandQueues, ctx->InputBuffer, CL_TRUE, 0, input_len*sizeof(uint8_t), input, 0, nullptr, nullptr)) != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clEnqueueWriteBuffer to fill input buffer.", err_to_str(ret));
        return OCL_ERR_API;
    }

    // Input
	const int cn0_kernel_offset = 0;//cn0KernelOffset(variant);
	
    if ((ret = OclLib::setKernelArg(ctx->Kernels[cn0_kernel_offset], 0, sizeof(cl_mem), &ctx->InputBuffer)) != CL_SUCCESS) {
        LOG_ERR(kSetKernelArgErr, err_to_str(ret), cn0_kernel_offset, 0);
        return OCL_ERR_API;
    }

    // output length
    if ((ret = OclLib::setKernelArg(ctx->Kernels[cn0_kernel_offset], 1, 8, &input_len)) != CL_SUCCESS) {
        LOG_ERR(kSetKernelArgErr, err_to_str(ret), cn0_kernel_offset,1);
        return OCL_ERR_API;
    }

	// output length
	if ((ret = OclLib::setKernelArg(ctx->Kernels[cn0_kernel_offset], 2, sizeof(cl_mem), &ctx->OutputBuffer)) != CL_SUCCESS) {
		LOG_ERR(kSetKernelArgErr, err_to_str(ret), cn0_kernel_offset, 2);
		return OCL_ERR_API;
	}

	
    return OCL_ERR_SUCCESS;
}

size_t XMRRunJob(GpuContext *ctx, cl_uint *HashOutput, xmrig::Variant variant)
{
    cl_int ret;
    cl_uint zero = 0;
	
    size_t g_intensity = ctx->rawIntensity;
    size_t w_size = OclCache::worksize(ctx, variant);
    // round up to next multiple of w_size
    size_t g_thd = ((g_intensity + w_size - 1u) / w_size) * w_size;
    // number of global threads must be a multiple of the work group size (w_size)
    assert(g_thd % w_size == 0);

    if ((ret = OclLib::enqueueWriteBuffer(ctx->CommandQueues, ctx->OutputBuffer, CL_FALSE,
		0, sizeof(cl_uint), &zero, 0, nullptr, nullptr)) != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clEnqueueWriteBuffer to fetch results.", err_to_str(ret));
        return OCL_ERR_API;
    }

    OclLib::finish(ctx->CommandQueues);

    size_t Nonce[2] = { ctx->Nonce, 1 }, gthreads[2] = { g_thd, 1 }, lthreads[2] = { w_size, 8 };
    const int cn0_kernel_offset = cn0KernelOffset(variant);

    if ((ret = OclLib::enqueueNDRangeKernel(ctx->CommandQueues, ctx->Kernels[0], 1, NULL, gthreads, lthreads, 0, nullptr, nullptr)) != CL_SUCCESS) {
        LOG_ERR("Error %s when calling clEnqueueNDRangeKernel for kernel %d.", err_to_str(ret), 0);
        return OCL_ERR_API;
    }

	if ((ret = OclLib::enqueueReadBuffer(ctx->CommandQueues, ctx->OutputBuffer, CL_TRUE,
		0, sizeof(cl_uint) *2, HashOutput, 0, nullptr, nullptr)) != CL_SUCCESS) {
		LOG_ERR("Error %s when calling clEnqueueWriteBuffer to fetch results.", err_to_str(ret));
		return OCL_ERR_API;
	}

	ctx->Nonce += (uint32_t)g_intensity;

    OclLib::finish(ctx->CommandQueues);
    // round up to next multiple of w_size


    return OCL_ERR_SUCCESS;
}


void ReleaseOpenCl(GpuContext* ctx)
{
    OclLib::releaseMemObject(ctx->InputBuffer);
    OclLib::releaseMemObject(ctx->OutputBuffer);


    OclLib::releaseProgram(ctx->Program);

    int kernel_count = sizeof(ctx->Kernels) / sizeof(ctx->Kernels[0]);
    for (int k = 0; k < kernel_count; ++k) {
        if (ctx->Kernels[k]) {
            OclLib::releaseKernel(ctx->Kernels[k]);
        }
    }

    OclLib::releaseCommandQueue(ctx->CommandQueues);
}


void ReleaseOpenClContext(cl_context opencl_ctx)
{
    OclLib::releaseContext(opencl_ctx);
}
