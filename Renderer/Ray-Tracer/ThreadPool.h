#pragma once
#include <thread>

#include <condition_variable>
#include <functional>
#include <vector>
#include <queue>
#include <iostream>

static class ThreadPool
{
public: 
	using Task = std::function<void()>;

	ThreadPool(short threads) {
		Start(threads);
	}

	~ThreadPool() {
		Stop();
	}
	
	void enqueue(Task task) {
		{
			// move is used to indicate that an object t may be "moved from",
			// i.e. allowing the efficient transfer of resources from t to another object. 
			std::unique_lock<std::mutex> lock(mMutex);
			mTasks.emplace(std::move(task));
		}

		// notify one waiting thread if there is one
		// can notify_all for all threads waiting
		mCondition.notify_one();
	}

	void updateThreads() {
		threads = std::thread::hardware_concurrency();
	}

	 short threads;
private:
	std::vector<std::thread> myThreads;
	std::vector<std::thread> mThreads;
	std::condition_variable mCondition;
	// mutex used by all threads - used to manage threads usage of data 
	std::mutex mMutex;
	std::queue<Task> mTasks;

	bool mRunning = false;

	void Stop() {
		{
			std::unique_lock<std::mutex> lock(mMutex);
			mRunning = true;
		}
		mCondition.notify_all();

		for (auto& thread : mThreads)
			if(thread.joinable())
				thread.join();
	};

	void Start(short threadCount) {
		for (short i = 0; i < threadCount; i++) {
			
			mThreads.emplace_back([=] {
					Task task;
				while (true) {
					{
						// lock thread so it cannot be used by other thread 
						// in order to prevent memory conflicts
						// unlocks after scope
						std::unique_lock<std::mutex> lock(mMutex);

						// wait for notification - spurious wake 
						// if queue empty put thread back to sleep else pop task
						mCondition.wait(lock, [=] {return mRunning || !mTasks.empty(); });

						if (mRunning && mTasks.empty())
							break;

						// return a referance of the next element in queue 
						task = std::move(mTasks.front());
						mTasks.pop();
					}
						task();
				}
				});
		}
	}
};

