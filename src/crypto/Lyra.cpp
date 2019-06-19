/* XMRig
 * Copyright 2010      Jeff Garzik <jgarzik@pobox.com>
 * Copyright 2012-2014 pooler      <pooler@litecoinpool.org>
 * Copyright 2014      Lucas Jones <https://github.com/lucasjones>
 * Copyright 2014-2016 Wolf9466    <https://github.com/OhGodAPet>
 * Copyright 2016      Jay D Dee   <jayddee246@gmail.com>
 * Copyright 2017-2019 XMR-Stak    <https://github.com/fireice-uk>, <https://github.com/psychocrypt>
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


#include <assert.h>
#include <sstream>


#include "common/cpu/Cpu.h"
#include "common/log/Log.h"
#include "common/net/Job.h"
#include "Mem.h"
#include "crypto/Lyra.h"
#include "crypto/lyra2.h"
#include "crypto/Lyra2_test.h"
#include "net/JobResult.h"


alignas(16) Lyra2_ctx *Lyra2::m_ctx = nullptr;
xmrig::Algo Lyra2::m_algorithm = xmrig::LYRA2;
xmrig::AlgoVerify Lyra2::m_av  = xmrig::VERIFY_HW_AES;


bool Lyra2::hash(const xmrig::Job &job, xmrig::JobResult &result, Lyra2_ctx *ctx)
{
    fn(job.algorithm().variant())(job.blob(), job.size(), result.result, &ctx, job.height());

    return *reinterpret_cast<uint64_t*>(result.result + 24) < job.target();
}

bool Lyra2::init(xmrig::Algo algorithm)
{
    m_algorithm = algorithm;
    m_av        = xmrig::Cpu::info()->hasAES() ? xmrig::VERIFY_HW_AES : xmrig::VERIFY_SOFT_AES;

    return selfTest();
}


template<xmrig::Algo ALGO, xmrig::Variant VARIANT>
static void cryptonight_single_hash_wrapper(const uint8_t *input, size_t size, uint8_t *output, Lyra2_ctx **ctx, uint64_t height)
{
    using namespace xmrig;
    lyra2_single_hash<ALGO, false, VARIANT>(input, size, output, ctx, height);

}


Lyra2::lyra2_hash_fun Lyra2::fn(xmrig::Algo algorithm, xmrig::AlgoVerify av, xmrig::Variant variant)
{
    using namespace xmrig;

    static const lyra2_hash_fun func_table[] = {
        lyra2_hash<xmrig::LYRA2, false, VARIANT_0>,
		lyra2_hash<xmrig::LYRA2, false, VARIANT_1>

    };

    return func_table[0];

}


bool Lyra2::selfTest() {
    using namespace xmrig;

    Mem::create(&m_ctx, m_algorithm, 1);

    if (m_algorithm == xmrig::LYRA2) {
        const bool rc = verify(VARIANT_0);

        if (!rc) {
            return rc;
        }
        return rc;
    }

    return false;
}


bool Lyra2::verify(xmrig::Variant variant)
{
    if (!m_ctx) {
        return false;
    }

    uint8_t output[32];

    lyra2_hash_fun func = fn(variant);
    if (!func) {
        return false;
    }

    func(test_input, 76, output, &m_ctx, 0);

    return memcmp(output, test_output_v1, 32) == 0;
}

