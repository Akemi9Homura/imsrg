#ifndef IMSRG_MEMORY_BUDGET_HH
#define IMSRG_MEMORY_BUDGET_HH

#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <limits>

namespace Commutator
{
namespace detail
{
// Budgets are in MiB. Invalid settings retain the default; zero disables a
// cache, or restricts a scratch-limited kernel to its required single worker.
inline size_t MemoryBudgetBytes(const char *name, size_t default_mib)
{
  const size_t mib=1024*1024;
  const char *value=std::getenv(name);
  if (!value || *value<'0' || *value>'9') return default_mib*mib;
  char *end=nullptr;
  errno=0;
  const auto parsed=std::strtoull(value,&end,10);
  if (errno || *end || parsed>std::numeric_limits<size_t>::max()/mib)
    return default_mib*mib;
  return static_cast<size_t>(parsed)*mib;
}
}
}
#endif
