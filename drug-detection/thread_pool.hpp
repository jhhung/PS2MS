#include <atomic>
#include <shared_mutex>
#include <future>
#include <chrono>
#include <queue>
#include <utility>
#include <functional>
#include <type_traits>

template <typename T>
struct pair_compare
{
    bool operator()(const T& lhs, const T& rhs)
    {
        return lhs.first < rhs.first;
    }
};

template <typename T>
class ThreadSafeQueue
{
  public:
    ThreadSafeQueue() {}

    void push(int priority, T&& value)
    {
        std::lock_guard<std::mutex> lock(mut);
        queue.emplace(priority, std::forward<T>(value));
        cond.notify_one();
    }

    T wait_and_pop()
    {
        std::unique_lock<std::mutex> lock(mut);
        cond.wait(lock, [this]{return !queue.empty();});
        T value = std::move(queue.top().second);
        queue.pop();
        return value;
    }

    bool try_pop(T& value)
    {
        std::lock_guard<std::mutex> lock(mut);
        if (queue.empty())
            return false;
        value = std::move(queue.top().second);
        queue.pop();
        return true;
    }

    unsigned int size()
    {
        return queue.size();
    }
  
  private:
    using ValueType = std::pair<int, T>;
    mutable std::mutex mut;
    std::priority_queue<ValueType, std::vector<ValueType>, pair_compare<ValueType>> queue;
    std::condition_variable cond;
};

template <typename T>
void empty_body(T result) {}

template<typename ReturnType, typename FuncType>
class ThreadPoolImpl
{
  public:
    ThreadPoolImpl(const std::function<FuncType>& func)
    : done          (false)
    , postprocess   (func)
    , worker_idle   (40, true)
    {
        try
        {
            const unsigned int thread_count = 40;//std::thread::hardware_concurrency();
            for (unsigned int i = 0; i < thread_count; ++i)
                threads.emplace_back(&ThreadPoolImpl::worker_thread, this, i);
        }
        catch(...)
        {
            done = true;
            throw;
        }
    }

    ~ThreadPoolImpl()
    {
        terminate_all_thread();
    }

    void submit(int priority, const std::function<ReturnType()>& func)
    {
        // while (task_queue.size() >= MAX_QUEUE_SIZE)
        //     std::this_thread::sleep_for(std::chrono::seconds(5));
        task_queue.push(priority, std::function<ReturnType()>(func));
    }

    unsigned int get_task_num()
    {
        return task_queue.size();
    }

    void terminate_all_thread()
    {
        done = true;
        for (auto& thread : threads)
        {
            if (thread.joinable())
                thread.join();
        }
    }

    bool is_idle()
    {
        bool r = true;
        for (const auto& item : worker_idle)
        {
            if (!item)
            {
                r = false;
                break;
            }
        }
        return r;
    }

  protected:
    void worker_thread(unsigned int num)
    {
        std::function<ReturnType()> task;
        while (!done)
        {
            if (task_queue.try_pop(task))
            {
                worker_idle[num] = false;
                if constexpr (std::is_void_v<ReturnType>)
                {
                    task();
                    postprocess();
                }
                else
                    postprocess(task());
            }
            else
            {
                worker_idle[num] = true;
                std::this_thread::yield();
            }
        }
    }

    std::atomic<bool> done;
    std::vector<bool> worker_idle;
    ThreadSafeQueue<std::function<ReturnType()>> task_queue;
    std::vector<std::thread> threads;
    std::function<FuncType> postprocess;
    const int MAX_QUEUE_SIZE = 2000000;
};

template <typename T>
class ThreadPool : public ThreadPoolImpl<T, void(T)>
{
  public:
    ThreadPool(std::function<void(T)> func = [](T) {})
    : ThreadPoolImpl<T, void(T)>(func)
    {}
};
template <>
class ThreadPool<void> : public ThreadPoolImpl<void, void()> 
{
  public:
    ThreadPool(std::function<void()> func = []() {})
    : ThreadPoolImpl<void, void()>(func)
    {}
};
