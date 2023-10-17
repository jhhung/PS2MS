#include <atomic>
#include <shared_mutex>
#include <future>
#include <chrono>
#include <queue>
#include <utility>
#include <functional>
#include <boost/shared_ptr.hpp>

template <typename T>
struct pair_compare
{
    bool operator()(const T& lhs, const T& rhs)
    {
        return lhs.first < rhs.first;
    }
};

class Result
{
  public:
    Result(std::string filename)
    : sum           (0)
    , num_predicted (0) 
    , outfile       (filename)
    {}

    void write_file(boost::shared_ptr<MolData> data)
    {
        std::lock_guard<std::mutex> lock(file_mut);
        data->writePredictedSpectraToMspFileStream( outfile );
    }

    void add_sum(const std::string& name, double score)
    {
        std::lock_guard<std::mutex> lock(mut);
        sum += score;
        ++num_predicted;
        std::cout << "Name: " << name << ", score: " << score << std::endl;
    }

    double get_avg_score()
    {
        std::lock_guard<std::mutex> lock(mut);
        return sum / num_predicted;
    }

    int size()
    {
        return num_predicted;
    }
  
  private:
    mutable std::mutex mut;
    mutable std::mutex file_mut;
    double sum;
    int num_predicted;
    std::ofstream outfile;
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

class ThreadPool
{
  public:
    ThreadPool()
    : done (false)
    {
        try
        {
            const unsigned int thread_count = std::thread::hardware_concurrency();
            for (unsigned int i = 0; i < thread_count; ++i)
                threads.emplace_back(&ThreadPool::worker_thread, this);
        }
        catch(...)
        {
            done = true;
            throw;
        }
    }

    ~ThreadPool()
    {
        terminate_all_thread();
    }

    template <typename T>
    void submit(int priority, T&& func)
    {
        task_queue.push(priority, std::forward<T>(func));
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

  private:
    void worker_thread()
    {
        std::function<void()> task;
        while (!done)
        {
            if (task_queue.try_pop(task))
                task();
            else
                std::this_thread::yield();
        }
    }

    std::atomic<bool> done;
    ThreadSafeQueue<std::function<void()>> task_queue;
    std::vector<std::thread> threads;
};